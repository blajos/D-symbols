var container, stats;
var camera, scene, renderer;
var splitting_system = new Array;
var points_labels = new Array;
var edges = new Array;
var g_transparency = true;
var g_opacity = 0.8;
var g_create_text = true;
var g_side = THREE.DoubleSide;
var g_update_controls = true;

init();
animate();

function new_simplex(z, num, maxnum, labels, colors) {
  // z=+1 vagy -1 (also vagy felso szimplex)
  // num: hanyadik (0-val kezdve, also es a felso kor kulon szamit)
  // maxnum: mennyibol (ez ugyanaz kell legyen alul-felul)
  // labels: Array of: { op: <elhagyott operacio>, simpleces: "<szokozzel
  //   felsorolva az erintett szimplexek rendezetten>" } a sorrend fontos: 1. az
  //   1-es operacio, 2. a 0-as operacio, 3-4. pedig felvaltva a 2-es, 3-as
  //   operacio.
  // parity: 0 vagy 1 az elozoekben emlitett felvaltva dolgozas elojelzese
  // colors: face colors: mindig a labels azonos indexuvel szemkozti szint
  // mondjuk meg.
  var parity = num % 2;
  var dszog = -2*Math.PI/maxnum;
  var szog = num*dszog;
  var simplex_v = [
    new THREE.Vector3(0,0,0),
	new THREE.Vector3(0,0,z), 
	new THREE.Vector3(Math.cos(szog),Math.sin(szog),0).multiplyScalar(Math.pow(Math.cos(dszog),parity)),
	new THREE.Vector3(Math.cos(szog+dszog),Math.sin(szog+dszog),0).multiplyScalar(Math.pow(Math.cos(dszog),1-parity))];

  if(edges.length == 0){
    var vertices=new Array;
    for (var z=-1;z<2;z+=2){
      var start=new THREE.Vector3(0,0,z);
      vertices.push(start);
      var dir1=new THREE.Vector3(Math.cos(szog),Math.sin(szog),0).multiplyScalar(parity*2);
      dir1.add(new THREE.Vector3(Math.cos(szog+dszog),Math.sin(szog+dszog),0).multiplyScalar((1-parity)*2));
      var dir2=new THREE.Vector3(Math.cos(szog+2*dszog),Math.sin(szog+2*dszog),0).multiplyScalar(parity*2);
      dir2.add(new THREE.Vector3(Math.cos(szog-dszog),Math.sin(szog-dszog),0).multiplyScalar((1-parity)*2));
      vertices.push(new THREE.Vector3().addVectors(start,dir1));
      vertices.push(new THREE.Vector3().addVectors(start,dir2));
      vertices.push(new THREE.Vector3().addVectors(start,new THREE.Vector3().addVectors(dir1,dir2)));
    }
    var geometry=new THREE.Geometry();
    geometry.vertices=vertices;
    var faces=[[0,1,3,2],[4,5,7,6],[0,4,5,1],[0,4,6,2],[3,7,5,1],[3,7,6,2]];
    for(var i=0;i<faces.length;i++){
      geometry.faces.push(new THREE.Face4(faces[i][0],faces[i][1],faces[i][2],faces[i][3]));
    }
    scene.add(new THREE.Mesh( 
	  geometry,
	  new THREE.MeshBasicMaterial( { color: 0x000000, side: g_side, wireframe: true } )));
  }

  for(var elhagy=0;elhagy < simplex_v.length; elhagy++){
    var simplex_g = new THREE.Geometry();

    for(var i=0;i<simplex_v.length; i++){
      if (i!=elhagy)
	simplex_g.vertices.push(simplex_v[i]);
    }
    for(var i=0;i<simplex_v.length; i++){
      for(var j=i+1;j<simplex_v.length; j++){
	if (i!=elhagy && j!=elhagy)
	  edges.push({from: simplex_v[i], to: simplex_v[j]});
      }
    }
    simplex_g.faces.push(new THREE.Face3(0,1,2));

    var simplex_m1 = new THREE.MeshBasicMaterial( {
      color: 0x000000,
	side: g_side,
	wireframe: true } );

    var simplex_m2 = new THREE.MeshBasicMaterial( {
      opacity: g_opacity,
	color: colors[elhagy],
	transparent: g_transparency,
	side: g_side,
	wireframe: false } );

    simplex_g.computeBoundingBox();
    simplex_g.computeFaceNormals();
    for (var j=0;j<3;j++){
      simplex_g.faces[0].vertexNormals.push(simplex_g.faces[0].normal.clone());
    }

    scene.add(new THREE.Mesh( simplex_g, simplex_m1 ));
    if (elhagy == 0) {
      scene.add(new THREE.Mesh( simplex_g, simplex_m2 ));
    }
  }


  //Text:
  for (var i = 0; i < simplex_v.length; i++){
    var j;
    for (j = 0; j < points_labels.length; j++){
      //console.log(points_labels[j], simplex_v[i], labels[i])
      if (points_labels[j].point.x == simplex_v[i].x &&
	  points_labels[j].point.y == simplex_v[i].y && 
	  points_labels[j].point.z == simplex_v[i].z && 
	  points_labels[j].op == labels[i].op &&
	  points_labels[j].simpleces == labels[i].simpleces){
	    //console.log("Helo :-)");
	    break;
	  }
    }
    if(j == points_labels.length){
      //console.log("Uj pont:", {point: simplex_v[i], op: labels[i].op,
      //simpleces: labels[i].simpleces});
      points_labels.push({point: simplex_v[i], op: labels[i].op, 
	simpleces: labels[i].simpleces});
      create_text(labels[i].op,labels[i].simpleces,simplex_v[i]);
    }
  }
}

function create_text(text, subscript, position){
  if (g_create_text != true) {
    return;
  }
  //console.log("Text");
  // Szovegek
  var textGeo = new THREE.TextGeometry( text, {

    size: 0.1,
      height: 0.001,
      curveSegments: 4,

      font: "helvetiker",
      weight: "normal",
      style: "normal",

      bevelThickness: 0.003,
      bevelSize: 0.002,
      bevelEnabled: true,

      material: 0,
      extrudeMaterial: 0.002

  });

  var tm = new THREE.MeshBasicMaterial( {
    color: 0x000000,
      opacity: 1,
      transparent: false,
      wireframe: false } );

  textGeo.computeBoundingBox();
  textGeo.computeVertexNormals();
  var textMesh = new THREE.Mesh( textGeo,
      tm );
  textMesh.position=position.clone();
  textMesh.position.z += 0.01;
  textMesh.rotation.x = Math.PI/2;

  scene.add(textMesh);

  var textGeo = new THREE.TextGeometry( subscript, {

    size: 0.04,
      height: 0.001,
      curveSegments: 4,

      font: "helvetiker",
      weight: "normal",
      style: "normal",

      bevelThickness: 0.002,
      bevelSize: 0.0008,
      bevelEnabled: true,

      material: 0,
      extrudeMaterial: 0.0008

  });

  textGeo.computeBoundingBox();
  textGeo.computeVertexNormals();
  var textMesh = new THREE.Mesh( textGeo,
      tm );
  textMesh.rotation.x = Math.PI/2;
  textMesh.rotation.z = -Math.PI/20;
  textMesh.position = position.clone();
  textMesh.position.x += 0.08;
  textMesh.position.z -= 0.01;

  scene.add(textMesh);
}

function create_splitting(labels){
  // labels: Array of: { op1: <elhagyott operacio>, simpleces1: "<szokozzel
  //   felsorolva az erintett szimplexek rendezetten>", op2: <tulso veg op>,
  //   simpleces2: <tulso veg simpleces> } 
  for (var i=0;i<splitting_system.length;i++)
    scene.remove(splitting_system[i]);
  splitting_system = new Array;

  for (var i = 0; i < labels.length; i++){
    // Find first and second points
    var p1=new Array;
    var p2=new Array;
    for (var j=0; j < points_labels.length; j++){
      if(labels[i].op1 == points_labels[j].op &&
	  labels[i].simpleces1 == points_labels[j].simpleces){
	    p1.push(points_labels[j].point);
	  }
      if(labels[i].op2 == points_labels[j].op &&
	  labels[i].simpleces2 == points_labels[j].simpleces){
	    p2.push(points_labels[j].point);
	  }
    }

    for (var i1=0; i1<p1.length; i1++){
      for (var i2=0; i2<p2.length; i2++){
	var is_really_edge=false;
	for (var j=0; j<edges.length; j++){
	  if((edges[j].from.equals(p1[i1]) && edges[j].to.equals(p2[i2])) ||
	     (edges[j].from.equals(p2[i2]) && edges[j].to.equals(p1[i1])))
	    is_really_edge=true;
	}
	if (is_really_edge){
	  var plane_g=new THREE.PlaneGeometry(0.1,0.1);
	  var plane_m1=new THREE.Mesh( plane_g, new THREE.MeshBasicMaterial( { color: 0xff0000}));
	  var plane_m2=new THREE.Mesh( plane_g, new THREE.MeshBasicMaterial( { color: 0x00ff00}));

	  var halfpos = new THREE.Vector3((p1[i1].x+p2[i2].x)/2,
	      (p1[i1].y+p2[i2].y)/2, (p1[i1].z+p2[i2].z)/2);
	  var normal = new THREE.Vector3().subVectors(p1[i1],halfpos).normalize();
	  var mx1 = new THREE.Matrix4(1,0,0,halfpos.x,
	                              0,1,0,halfpos.y,
				      0,0,normal.z,halfpos.z,
				      0,0,0,1);
          if (normal.x != 0 || normal.y != 0){
	    var b1=new THREE.Vector3().crossVectors(normal, new
		THREE.Vector3(0,0,1));
	    var b2=new THREE.Vector3().crossVectors(normal, b1);
	    mx1.elements[0]=b1.x;
	    mx1.elements[4]=b1.y;
	    mx1.elements[8]=b1.z;
	    mx1.elements[1]=b2.x;
	    mx1.elements[5]=b2.y;
	    mx1.elements[9]=b2.z;
	    mx1.elements[2]=normal.x;
	    mx1.elements[6]=normal.y;
	    mx1.elements[10]=normal.z;
	  }
	  mx1.transpose();
          plane_m1.applyMatrix(mx1);
	  plane_m1.position=halfpos;
	  mx1.elements[8]=-mx1.elements[8];
	  mx1.elements[9]=-mx1.elements[9];
	  mx1.elements[10]=-mx1.elements[10];
          plane_m2.applyMatrix(mx1);
	  plane_m2.position=halfpos;
	  splitting_system.push(plane_m1);
	  scene.add(plane_m1);
	  splitting_system.push(plane_m2);
	  scene.add(plane_m2);
	}
      }
    }
  }
}

function init() {

  container = document.getElementById( 'container' );

  //stats = new Stats();

  //camera = new THREE.OrthographicCamera(-3, 3, -2, 2, 0.1, 10);
  camera = new THREE.PerspectiveCamera( 40, window.innerWidth /
      window.innerHeight, 0.0001, 1000 );
  camera.position.y = -4;
  camera.position.z = 2;
  camera.up.x = 0;
  camera.up.y = 0;
  camera.up.z = 1;
  controls = new THREE.TrackballControls( camera );

  scene = new THREE.Scene();

  if ( ! Detector.webgl ){
    Detector.addGetWebGLMessage();
    renderer = new THREE.SVGRenderer();
    g_transparency = false;
    g_opacity = 1;
    g_create_text = false;
    //g_side = THREE.FrontSide;
    g_update_controls = false;
  }
  else {
    renderer = new THREE.WebGLRenderer({antialias: true});
  }
  renderer.setSize( window.innerWidth, window.innerHeight );

  container.appendChild( renderer.domElement );

  //stats = new Stats();
  //stats.domElement.style.position = 'absolute';
  //stats.domElement.style.top = '0px';
  //container.appendChild( stats.domElement );

  window.addEventListener( 'resize', onWindowResize, false );

  controls.update();
}

function onWindowResize() {

  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();

  renderer.setSize( window.innerWidth/2, window.innerHeight/2 );
}

function animate() {

  requestAnimationFrame( animate );

  render();
  //stats.update();

}

function render() {

  if (g_update_controls == true){
    controls.update();
  }
  renderer.render( scene, camera );

}


