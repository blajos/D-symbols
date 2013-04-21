var container, stats;
var camera, scene, renderer;
var splitting_system;
var points_labels = new Array;

init();
animate();

function new_simplex(z, num, maxnum, parity, labels) {
  // z=+1 vagy -1 (also vagy felso szimplex)
  // num: hanyadik (0-val kezdve, also es a felso kor kulon szamit)
  // maxnum: mennyibol (ez ugyanaz kell legyen alul-felul)
  // labels: Array of: { op: <elhagyott operacio>, simpleces: "<szokozzel
  //   felsorolva az erintett szimplexek rendezetten>" } a sorrend fontos: 1. az
  //   1-es operacio, 2. a 0-as operacio, 3-4. pedig felvaltva a 2-es, 3-as
  //   operacio.
  // parity: 0 vagy 1 az elozoekben emlitett felvaltva dolgozas elojelzese
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
  for (var i = 0; i < simplex_v.length; i++){
    points_labels.push({point: simplex_v[i], op: labels[i].op, 
      simpleces: labels[i].simpleces});
  }
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
    create_text(labels[i].op,labels[i].simpleces,simplex_v[i]);
  }
}

function create_text(text, subscript, position){
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
  textMesh.position=position;
  textMesh.position.z += 0.01;
  textMesh.rotation.x = Math.PI/2;

  scene.add(textMesh);

  var textGeo = new THREE.TextGeometry( subscript, {

    size: 0.02,
      height: 0.001,
      curveSegments: 4,

      font: "helvetiker",
      weight: "normal",
      style: "normal",

      bevelThickness: 0.001,
      bevelSize: 0.0004,
      bevelEnabled: true,

      material: 0,
      extrudeMaterial: 0.0004

  });

  textGeo.computeBoundingBox();
  textGeo.computeVertexNormals();
  var textMesh = new THREE.Mesh( textGeo,
      tm );
  textMesh.rotation.x = Math.PI/2;
  textMesh.position = position;
  textMesh.position.x += 0.08;
  textMesh.position.z += 0.01;

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
    var p1=null;
    var p2=null;
    for (var j=0; j < points_labels.length; j++){
      if(p1 == null && labels[i].op1 == points_labels[j].op &&
	  labels[i].simpleces1 == points_labels[j].simpleces){
	    p1 = points_labels[j].point;
	  }
      if(p2 == null && labels[i].op2 == points_labels[j].op &&
	  labels[i].simpleces2 == points_labels[j].simpleces){
	    p2 = points_labels[j].point;
	  }
    }

    var moveto = [0,0,0];
    for (var k=0; k<p1.length; k++){
      moveto[k]=(p1[k]+p2[k])/2;
    }

    var vertex = new THREE.Vector3();
    vertex.x = moveto[0];
    vertex.y = moveto[1];
    vertex.z = moveto[2];
    geometry.vertices.push( vertex );
  }
  splitting_system = new THREE.ParticleSystem(geometry, new
      THREE.ParticleBasicMaterial({color: 0x0f0f0f, size: 0.2, opacity: 0.8}) );
  scene.add(splitting_system);
}

function init() {
  if ( ! Detector.webgl ) Detector.addGetWebGLMessage();

  container = document.getElementById( 'container' );

  stats = new Stats();

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
  renderer.setSize( window.innerWidth, window.innerHeight );

  container.appendChild( renderer.domElement );

  stats = new Stats();
  stats.domElement.style.position = 'absolute';
  stats.domElement.style.top = '0px';
  container.appendChild( stats.domElement );

  window.addEventListener( 'resize', onWindowResize, false );

}

function onWindowResize() {

  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();

  renderer.setSize( window.innerWidth, window.innerHeight );
}

function animate() {

  requestAnimationFrame( animate );

  render();
  stats.update();

}

function render() {

  controls.update();
  renderer.render( scene, camera );

}


