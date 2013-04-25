var container, stats;
var camera, scene, renderer;
var splitting_system;
var points_labels = new Array;

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
  var dszog = 2*Math.PI/maxnum;
  var szog = num*dszog;
  var simplex_v = [
    new THREE.Vector3(0,0,0),
	new THREE.Vector3(0,0,z), 
	new THREE.Vector3(Math.pow(Math.cos(dszog),parity)*Math.cos(szog), 
	    Math.pow(Math.cos(dszog),parity)*Math.sin(szog),
	    0),
	new THREE.Vector3(Math.pow(Math.cos(dszog),1-parity)*Math.cos(szog+dszog),
	    Math.pow(Math.cos(dszog),1-parity)*Math.sin(szog+dszog),
	    0)];
  var simplex_faces = [[1,2,3],
      [0,2,3],
      [0,1,3],
      [0,1,2]];

  var simplex_g = new THREE.Geometry();
  var simplex_materials = new Array;

  for(var i=0;i<simplex_v.length; i++){
    simplex_g.vertices.push(simplex_v[i]);

    simplex_materials.push(new THREE.MeshBasicMaterial( {
      //opacity: 0.5,
      color: colors[i],
      //transparent: true,
      //vertexColors: THREE.FaceColors,
      side: THREE.DoubleSide,
      shading: THREE.FlatShading, 
      overdraw: true,
      wireframe: false } ));

    //console.log("Color:",colors[i]);
    var currface=new THREE.Face3(simplex_faces[i][0],simplex_faces[i][1],simplex_faces[i][2]);
    currface.materialIndex = i;
    currface.color.setHex(colors[i]);
    simplex_g.faces.push(currface);
    var faceuv = [
                    new THREE.UV(0,1),
                    new THREE.UV(1,1),
                    new THREE.UV(1,0)
                ];
    simplex_g.faceUvs[0].push(new THREE.UV(0,1));
    simplex_g.faceVertexUvs[0].push(faceuv);
  }
  //simplex_g.materials = simplex_materials;

  var simplex_m1 = new THREE.MeshBasicMaterial( {
    color: 0x000000,
      wireframe: true } );

  var simplex_m2 = new THREE.MeshBasicMaterial( {
    opacity: 0.8,
    color: colors[0],
    transparent: true,
    //vertexColors: THREE.FaceColors,
    side: THREE.DoubleSide,
    //shading: THREE.FlatShading, 
    //overdraw: true,
    wireframe: false } );

  simplex_g.computeBoundingBox();
  //simplex_g.computeCentroids();
  simplex_g.computeFaceNormals();
  for(var i=0;i<simplex_g.faces.length;i++){
    for (var j=0;j<3;j++){
      simplex_g.faces[i].vertexNormals.push(simplex_g.faces[i].normal.clone());
    }
  }
  //simplex_g.computeMorphNormals();
  //simplex_g.computeTangents();
  //simplex_g.computeLineDistances();


  scene.add(new THREE.Mesh( simplex_g, simplex_m1 ));
  console.log("Geometry:",simplex_g);
  //scene.add(new THREE.Mesh( simplex_g, new THREE.MeshFaceMaterial(simplex_g.materials) ));
  scene.add(new THREE.Mesh( simplex_g, simplex_m2 ));

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
  try {
    scene.remove(splitting_system);
  }
  catch (err){}

  var geometry = new THREE.Geometry();
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
	var vertex = new THREE.Vector3();
	vertex.x = (p1[i1].x+p2[i2].x)/2;
	vertex.y = (p1[i1].y+p2[i2].y)/2;
	vertex.z = (p1[i1].z+p2[i2].z)/2;
	geometry.vertices.push( vertex );
      }
    }
  }
  splitting_system = new THREE.ParticleSystem(geometry, new
      THREE.ParticleBasicMaterial({color: 0x0f0f0f, size: 0.2, opacity: 0.8}) );
  scene.add(splitting_system);
}

function init() {
  if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

  container = document.getElementById( 'container' );

  //stats = new Stats();

  camera = new THREE.PerspectiveCamera( 70,
      window.innerWidth / window.innerHeight, 0.1, 10 );
  camera.position.y = -2;
  camera.position.z = 1;
  camera.up.x = 0;
  camera.up.y = 0;
  camera.up.z = 1;
  controls = new THREE.TrackballControls( camera );

  scene = new THREE.Scene();

  renderer = new THREE.WebGLRenderer({antialias: true});
  //renderer = new THREE.CanvasRenderer();
  renderer.setSize( window.innerWidth/2, window.innerHeight/2 );

  container.appendChild( renderer.domElement );

  //stats = new Stats();
  //stats.domElement.style.position = 'absolute';
  //stats.domElement.style.top = '0px';
  //container.appendChild( stats.domElement );

  window.addEventListener( 'resize', onWindowResize, false );

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

  controls.update();
  renderer.render( scene, camera );

}


