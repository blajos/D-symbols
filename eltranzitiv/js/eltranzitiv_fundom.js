var container, stats;
var camera, scene, renderer;
var splitting_system;
var points_labels = new Array;

init();
animate();

function new_simplex(z, num, maxnum, labels) {
  // z=+1 vagy -1 (also vagy felso szimplex)
  // num: hanyadik (0-val kezdve, also es a felso kor kulon szamit)
  // maxnum: mennyibol (ez ugyanaz kell legyen alul-felul)
  // labels: Array of: { op: <elhagyott operacio>, simpleces: "<szokozzel
  //   felsorolva az erintett szimplexek rendezetten>" } a sorrend fontos: 1. az
  //   1-es operacio, 2. a 0-as operacio, 3-4. pedig felvaltva a 2-es, 3-as
  //   operacio.
  // parity: 0 vagy 1 az elozoekben emlitett felvaltva dolgozas elojelzese
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
  var simplex_faces = [[2,1,0],[0,3,2],[1,3,0],[2,3,1]];
//  var simplex_g = new THREE.PolyhedronGeometry(
//      simplex_v, simplex_faces);
  var simplex_g = new THREE.ConvexGeometry( simplex_v );
  var simplex_c = Math.random() * 0xaaaaaa;
  var simplex_m1 = new THREE.MeshBasicMaterial( {
    color: 0x000000,
      wireframe: true } );

  var simplex_m2 = new THREE.MeshBasicMaterial( {
    color: simplex_c,
      opacity: 0.5,
      transparent: true,
      wireframe: false } );

  var simplex = new THREE.Mesh( simplex_g,
      simplex_m1 );
  scene.add(simplex);
  scene.add(new THREE.Mesh( simplex_g,
	simplex_m2 ));

  //Text:
  for (var i = 0; i < simplex_v.length; i++){
    var j;
    for (j = 0; j < points_labels.length; j++){
      if (points_labels[j].point == simplex_v[i] &&
	  points_labels[j].op == labels[i].op &&
	  points_labels[j].simpleces == labels[i].simpleces){
	    break;
	  }
    }
    if(j == points_labels.length){
    points_labels.push({point: simplex_v[i], op: labels[i].op, 
      simpleces: labels[i].simpleces});
    create_text(labels[i].op,labels[i].simpleces,simplex_v[i]);
    }
  }
}

function create_text(text, subscript, position){
  var canvas1 = document.createElement('canvas');
  var context1 = canvas1.getContext('2d');
  context1.font = "Bold 40px Arial";
  context1.fillStyle = "rgba(15,15,15,1)";
  context1.fillText(text, 0, 50);

  // canvas contents will be used for a texture
  var texture1 = new THREE.Texture(canvas1) 
    texture1.needsUpdate = true;

  var material1 = new THREE.MeshBasicMaterial( {map: texture1, side:THREE.DoubleSide } );
  material1.transparent = true;

  var mesh1 = new THREE.Mesh(
      new THREE.PlaneGeometry(canvas1.width, canvas1.height),
      material1
      );
  mesh1.scale.x = mesh1.scale.y = mesh1.scale.z = 0.001;
  mesh1.position=position.clone();
  mesh1.position.z += 0.01;
  mesh1.rotation.x = Math.PI/2;
  scene.add( mesh1 );

  var canvas2 = document.createElement('canvas');
  var context2 = canvas2.getContext('2d');
  context2.font = "Bold 16px Arial";
  context2.fillStyle = "rgba(15,15,15,1)";
  context2.fillText(subscript, 0, 50);

  // canvas contents will be used for a texture
  var texture2 = new THREE.Texture(canvas2) 
    texture2.needsUpdate = true;

  var material1 = new THREE.MeshBasicMaterial( {map: texture2, side:THREE.DoubleSide } );
  material1.transparent = true;

  var mesh2 = new THREE.Mesh(
      new THREE.PlaneGeometry(canvas2.width, canvas2.height),
      material1
      );
  console.log("Canvas2:", canvas2.width, canvas2.height);
  mesh2.scale.x = mesh2.scale.y = mesh2.scale.z = 0.001;
  mesh2.position=position.clone();
  mesh2.position.x += 0.08;
  mesh2.position.z -= 0.01;
  mesh2.rotation.x = Math.PI/2;
  scene.add( mesh2 );
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


