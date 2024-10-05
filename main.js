// Imports
import * as THREE from 'https://cdn.skypack.dev/three@0.128.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.skypack.dev/three@0.128.0/examples/jsm/controls/OrbitControls.js';

// Constants
const DEG_TO_RAD = Math.PI / 180;
const AU_SCALE_FACTOR = 50;
const TA_TIME_SCALE_FACTOR = 0.0001; // This will not be needed when the true anomaly code is included

const DEFAULT_MESH_N = 32;
const ORBIT_MESH_POINTS = 128;
const ORBIT_SEGMENT_CONST = 2 * Math.PI / ORBIT_MESH_POINTS;

const NEO_ORBIT_COLOR = 0x1e90FF;
const NEO_COLOR = 0xFFFFFF;
const NEO_RADIUS = 0.5;
const MAX_VISIBLE_NEOS = 50;

const MOUSE_MIN_MOVE_CLICK = 0.01;

// FPS control
const targetFPS = 60; // Target frames per second
const frameInterval = 1000 / targetFPS; // Time per frame in milliseconds
let lastFrameTime = 0; // Tracks the last frame's timestamp

// Setup Scene
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 10000);
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

camera.position.set(100, 100, 100);
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

function createOrbit(orbitParams, color) {
    const cosNode = Math.cos(orbitParams.node);
    const sinNode = Math.sin(orbitParams.node);
    const cosPeri = Math.cos(orbitParams.peri);
    const sinPeri = Math.sin(orbitParams.peri);
    const cosInc = Math.cos(orbitParams.inc);
    const sinInc = Math.sin(orbitParams.inc);

    const row1 = [cosPeri * cosNode - cosInc * sinPeri * sinNode, -cosNode * sinPeri - cosInc * cosPeri * sinNode, sinInc * sinNode];
    const row2 = [cosPeri * sinNode + cosInc * cosNode * sinPeri, -sinPeri * sinNode + cosInc * cosPeri * cosNode, -sinInc * cosNode];
    const row3 = [sinInc * sinPeri, sinInc * cosPeri, cosInc];
    const matrix = [row1, row2, row3];

    orbitParams['transformMatrix'] = matrix;

    const points = [];
    const b = orbitParams.a * Math.sqrt(1 - orbitParams.e ** 2); // Semi-minor axis

    for (let i = 0; i <= ORBIT_MESH_POINTS; i++) {
        const eccentric_anomaly = ORBIT_SEGMENT_CONST * i; // Angle
        const xOrb = orbitParams.a * (Math.cos(eccentric_anomaly) - orbitParams.e);
        const yOrb = b * Math.sin(eccentric_anomaly);

        const xCamera = matrix[0][0] * xOrb + matrix[0][1] * yOrb;
        const yCamera = matrix[1][0] * xOrb + matrix[1][1] * yOrb;
        const zCamera = matrix[2][0] * xOrb + matrix[2][1] * yOrb;

        points.push(new THREE.Vector3(xCamera, zCamera, -yCamera));
    }

    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    const material = new THREE.LineBasicMaterial({ color: color });
    return new THREE.Line(geometry, material);
}

function getOrbitPosition(a, e, trueAnomaly, matrix) {
    const cosTA = Math.cos(trueAnomaly);
    const sinTA = Math.sin(trueAnomaly);
    const radius = a * (1 - e * e) / (1 + e * cosTA);

    const xOrb = radius * cosTA;
    const yOrb = radius * sinTA;

    const xCamera = matrix[0][0] * xOrb + matrix[0][1] * yOrb;
    const yCamera = matrix[1][0] * xOrb + matrix[1][1] * yOrb;
    const zCamera = matrix[2][0] * xOrb + matrix[2][1] * yOrb;

    return new THREE.Vector3(xCamera, zCamera, -yCamera);
}

function addSun() {
    // add Sun Texture
    const sunTextureLoader = new THREE.TextureLoader();
    const sunTexture = sunTextureLoader.load(
        'assets/body_textures/8k_sun.jpg'
    );

    const geometry = new THREE.SphereGeometry(1, DEFAULT_MESH_N, DEFAULT_MESH_N);
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
        orbitParams.a *= AU_SCALE_FACTOR;
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

        const orbit = createOrbit(orbitParams, planetData.renderParams.color);
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
        orbitParams.a *= AU_SCALE_FACTOR;

        const geometry = new THREE.SphereGeometry(NEO_RADIUS, DEFAULT_MESH_N / 2, DEFAULT_MESH_N / 2);
        const material = new THREE.MeshBasicMaterial({ color: NEO_COLOR });
        const neoMesh = new THREE.Mesh(geometry, material);
        neoMeshes[neoName] = neoMesh;

        const orbit = createOrbit(orbitParams, NEO_ORBIT_COLOR);
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