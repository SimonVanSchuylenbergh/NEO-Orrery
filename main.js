// Imports
import * as THREE from 'https://cdn.skypack.dev/three@0.128.0/build/three.module.js';
import { GLTFLoader } from 'https://cdn.jsdelivr.net/npm/three@0.114/examples/jsm/loaders/GLTFLoader.js';
import { OrbitControls } from 'https://cdn.skypack.dev/three@0.128.0/examples/jsm/controls/OrbitControls.js';
import { createOrbit, getOrbitPosition } from './orbits.js'

// Constants
const DEG_TO_RAD = Math.PI / 180;
const TA_TIME_SCALE_FACTOR = 0.0001; // This will not be needed when the true anomaly code is included

const DEFAULT_MESH_N = 32;
const ORBIT_MESH_POINTS = 128;

const NEO_ORBIT_COLOR = 0x1e90FF;
const NEO_COLOR = 0xFFFFFF;
const NEO_RADIUS = 0.01;
const MAX_VISIBLE_NEOS = 50;

const MOUSE_MIN_MOVE_CLICK = 0.01;

// FPS control
const targetFPS = 60; // Target frames per second
const frameInterval = 1000 / targetFPS; // Time per frame in milliseconds
let lastFrameTime = 0; // Tracks the last frame's timestamp

// Setup Scene
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer({antialias:true, canvas: document.getElementById("orreryCanvas")});
renderer.setSize(window.innerWidth, window.innerHeight);

// Load skybox texture
scene.background = new THREE.CubeTextureLoader().load([
    'assets/px.png', // Right
    'assets/nx.png', // Left
    'assets/py.png', // Top
    'assets/ny.png', // Bottom
    'assets/pz.png', // Front
    'assets/nz.png'  // Back
]);

// Set camera position
camera.position.z = 1;

camera.position.set(2, 2, 2);
camera.lookAt(0, 0, 0);

// Add lighting
const light = new THREE.AmbientLight(0x404040, 0.5); // Soft white light
scene.add(light);

// Setup Controls
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.update();

// Setup Event Listeners
window.addEventListener("resize", () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
});

const mouseDownXY = new THREE.Vector2();
const mouseUpXY = new THREE.Vector2();
const raycaster = new THREE.Raycaster();

document.addEventListener('mousedown', (event) => {
    mouseDownXY.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouseDownXY.y = -(event.clientY / window.innerHeight) * 2 + 1;

    console.log(mouseUpXY.x, mouseUpXY.y, "DOWN");
});

document.addEventListener('mouseup', (event) => {
    mouseUpXY.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouseUpXY.y = -(event.clientY / window.innerHeight) * 2 + 1;

    console.log(mouseUpXY.x, mouseUpXY.y, "UP");

    //if (mouseDownXY.sub(mouseUpXY)) {

    //}

    //update the picking ray with the camera and pointer position
    raycaster.setFromCamera(mouseUpXY, camera);

    // alculate objects intersecting the picking ray
    const intersects = raycaster.intersectObjects(scene.children);

    for (let i = 0; i < intersects.length; i ++ ) {
        intersects[i].object.material.color.set(0x00ff00); 
    }
});


document.addEventListener("click", (event) => {
    mouseDownXY.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouseDownXY.y = -(event.clientY / window.innerHeight) * 2 + 1;

    console.log(mouseDownXY.x, mouseDownXY.y, "CLICK");
});


// Functions
async function readJSON(filePath) {
    try {
        const response = await fetch(filePath);
        if (!response.ok) { throw new Error('Network response was not ok'); }
        const data = await response.json(); // Parse the JSON from the response
        return data; // Return the parsed data
    } catch (error) { console.error('There was a problem trying to read ' + filePath + ':', error); }
}

function addSun() {
    // add Sun Texture
    const sunTextureLoader = new THREE.TextureLoader();
    const sunTexture = sunTextureLoader.load(
        'assets/body_textures/8k_sun.jpg'
    );

    const geometry = new THREE.SphereGeometry(0.02, DEFAULT_MESH_N, DEFAULT_MESH_N);
    const material = new THREE.MeshBasicMaterial({map: sunTexture});
    const sunMesh = new THREE.Mesh(geometry, material);
    scene.add(sunMesh);
}

function initializePlanets() {
    for (const [planetName, planetData] of Object.entries(planets)) {
        const orbitParams = planetData.orbitParams;
        orbitParams.inc *= DEG_TO_RAD;
        orbitParams.node *= DEG_TO_RAD;
        orbitParams.peri *= DEG_TO_RAD;
        orbitParams.ma *= DEG_TO_RAD;
        // get planet texture
        const planetTextureName = planetData.renderParams.texture;
        const planetTextureLoader = new THREE.TextureLoader();
        const planetTexture = planetTextureLoader.load(
            'assets/body_textures/' + planetTextureName
        );

        const geometry = new THREE.SphereGeometry(planetData.renderParams.radius, DEFAULT_MESH_N, DEFAULT_MESH_N);
        const material = new THREE.MeshBasicMaterial({map: planetTexture}); // add texture
        const planetMesh = new THREE.Mesh(geometry, material);
        planetMeshes[planetName] = planetMesh;

        const orbit = createOrbit(orbitParams, planetData.renderParams.color, ORBIT_MESH_POINTS);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, 0, orbitParams.transformMatrix);
        planetMesh.position.set(pos.x, pos.y, pos.z);

        scene.add(orbit);
        scene.add(planetMesh);
    }
}

function initializeNeos() {
    let i = 0;
    for (const [neoName, neoData] of Object.entries(neos)) {
        const orbitParams = neoData.orbitParams;
        orbitParams.inc *= DEG_TO_RAD;
        orbitParams.node *= DEG_TO_RAD;
        orbitParams.peri *= DEG_TO_RAD;
        orbitParams.ma *= DEG_TO_RAD;

        const geometry = new THREE.SphereGeometry(NEO_RADIUS, DEFAULT_MESH_N / 2, DEFAULT_MESH_N / 2);
        const material = new THREE.MeshBasicMaterial({ color: NEO_COLOR });
        const neoMesh = new THREE.Mesh(geometry, material);
        neoMeshes[neoName] = neoMesh;

        const orbit = createOrbit(orbitParams, NEO_ORBIT_COLOR, ORBIT_MESH_POINTS);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, 0, orbitParams.transformMatrix);
        neoMesh.position.set(pos.x, pos.y, pos.z);

        scene.add(orbit);
        scene.add(neoMesh);

        i += 1;
        if (i == MAX_VISIBLE_NEOS) { break };
    }
}

// Data
let sunMesh;
const planetMeshes = {};
const neoMeshes = {};

const planets = await readJSON('data/planet_data.json');
const neos = await readJSON('data/risk_list_neo_data.json');

addSun();
initializePlanets(); // Initialize planets once
initializeNeos(); // Initialize NEOs once
// console.log(planets.Saturn);

// Add Saturn's Rings:
// const saturnData = planets.Saturn
// const ringLoader = new GLTFLoader();
// ringLoader.load('assets/body_textures/Saturn_rings.glb', function (gltf) {
//     const ringTexture = gltf.scene.children[0].material.map; // Extract texture from .glb
//     // Load Saturn again to put its rings on it
//     const ringGeometry = new THREE.SphereGeometry(0.02, DEFAULT_MESH_N, DEFAULT_MESH_N);  // Create a sphere
//     const ringMaterial = new THREE.MeshBasicMaterial({ map: ringTexture });  // Apply the loaded texture
//     const ringSphere = new THREE.Mesh(ringGeometry, ringMaterial);  // Create a mesh with the geometry and material
//     // Orbits
//     // const ringOrbit = createOrbit(saturnData.orbitParams, saturnData.renderParams.color, ORBIT_MESH_POINTS);
//     // const ringPos = getOrbitPosition(saturnData.orbitParams.a, saturnData.orbitParams.e, 0, saturnData.orbitParams.transformMatrix);
//     ringSphere.position.set(0, 0, 0);
//     // add to scene
//      //scene.add(ringOrbit);
//     scene.add(ringSphere);
// }, undefined, function (error) {
//     console.error('An error occurred loading the GLB:', error);
// });

// Animation loop with FPS control
function animate(time) {
    requestAnimationFrame(animate);

    // Limit frame rate
    const deltaTime = time - lastFrameTime;
    if (deltaTime < frameInterval) {
        return; // Skip frame if too soon
    }
    lastFrameTime = time;

    let currentTime = Date.now() * TA_TIME_SCALE_FACTOR;

    // Update planet positions
    for (const [planetName, neoData] of Object.entries(planets)) {
        const orbitParams = neoData.orbitParams;
        const trueAnomaly = currentTime;
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, trueAnomaly, orbitParams.transformMatrix);
        planetMeshes[planetName].position.set(pos.x, pos.y, pos.z);
    }

    // Update NEO positions
    let i = 0;
    for (const [neoName, neoData] of Object.entries(neos)) {
        const orbitParams = neoData.orbitParams;
        const trueAnomaly = currentTime;
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, trueAnomaly, orbitParams.transformMatrix);
        neoMeshes[neoName].position.set(pos.x, pos.y, pos.z);

        i += 1;
        if (i == MAX_VISIBLE_NEOS) { break };
    }

    controls.update();
    renderer.render(scene, camera);
}

// Start animation loop
requestAnimationFrame(animate);