// Imports
import * as THREE from 'https://cdn.skypack.dev/three@0.124.0/build/three.module.js';
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

const MOUSE_MIN_MOVE_CLICK = 0.005;

// FPS control
const targetFPS = 60; // Target frames per second
const frameInterval = 1000 / targetFPS; // Time per frame in milliseconds
let lastFrameTime = 0; // Tracks the last frame's timestamp

// Setup Scene
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.01, 1000);
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

const mouseDownXY = new THREE.Vector2(-10, -10);
const mouseUpXY = new THREE.Vector2(-10, -10);
const raycaster = new THREE.Raycaster(); //ray through the screen at the location of the mouse pointer (when the mouse is released)
raycaster.params.Line = {threshold: 0.1}; //needs to depend on zoom level
raycaster.params.Points = {threshold: 0.1}
let highlightedObj = null;
let prevColor = 0;
let moved = false;
let stackedObjIndex = 0; //the index of the array of all the uniquely selected objects that the casted ray intersected

//activates when the mouse is pressed down
document.addEventListener('pointerdown', (event) => {
    mouseDownXY.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouseDownXY.y = -(event.clientY / window.innerHeight) * 2 + 1;

    if (Math.abs(mouseDownXY.x - mouseUpXY.x) < MOUSE_MIN_MOVE_CLICK && Math.abs(mouseDownXY.y - mouseUpXY.y) < MOUSE_MIN_MOVE_CLICK) {
        moved = false; //mouse movement was small enough to not count as a move
    }
    else {
        moved = true;
        stackedObjIndex = 0;
    }
});

//activates when the mouse is released
document.addEventListener('pointerup', (event) => {
    mouseUpXY.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouseUpXY.y = -(event.clientY / window.innerHeight) * 2 + 1;

    if (Math.abs(mouseDownXY.x - mouseUpXY.x) < MOUSE_MIN_MOVE_CLICK && Math.abs(mouseDownXY.y - mouseUpXY.y) < MOUSE_MIN_MOVE_CLICK) { //didn't move mouse
        //update the picking ray with the camera and pointer position
        raycaster.setFromCamera(mouseUpXY, camera);

        // calculate objects intersecting the picking ray
        const allSelectedObjs = raycaster.intersectObjects(scene.children);

        //remove duplicate intersections of the same orbit
        const seen = new Set();
        const selectedObjs = allSelectedObjs.filter(item => {
            if (!seen.has(item.object.uuid)) {
                seen.add(item.object.uuid);
                return true; // unique uuid
            }
            return false; // same uuid
        });

        if (highlightedObj != null){ //clicking on the background deselects the current object (if there is one)
            highlightedObj.material.color.set(prevColor);
            highlightedObj = null;
        }
        
        if (selectedObjs.length != 0){
            if (!moved) { stackedObjIndex = (stackedObjIndex + 1) % selectedObjs.length; }
            
            if (selectedObjs[stackedObjIndex].object.type == 'Line') { 
                highlightedObj = selectedObjs[stackedObjIndex].object; //save the highlighted object
                prevColor = highlightedObj.material.color.getHex(); //save the highlighted object's previous color
                highlightedObj.material.color.set(0x00ff00); //highlight the select object if it is an orbit
            }
        }
    }
    else { //moved mouse
        moved = true;
        stackedObjIndex = 0;
    }
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
        // check if planet is saturn's rings
        // if so, make it a ring geometry with specified parameters -- otherwise, make it a spher egeometry
        // console.log(planetName)
        if (planetName == 'rings'){
            const geometry = new THREE.RingGeometry(planetData.renderParams.innerRadius, 
                planetData.renderParams.outerRadius, 64);
            const material = new THREE.MeshBasicMaterial({
                map: planetTexture,
                side: THREE.DoubleSide,
                transparent: true,
                opacity: 0.7
            });

            // Create the mesh
            var mesh = new THREE.Mesh(geometry, material);
            mesh.rotation.x = Math.PI / 2; //Rotate the rings to be flat
        }
        // if not ring do sphere
        else {
            const geometry = new THREE.SphereGeometry(planetData.renderParams.radius, DEFAULT_MESH_N, DEFAULT_MESH_N);
            const material = new THREE.MeshBasicMaterial({map: planetTexture}); // add texture
            // console.log(planetName)
            // Create the mesh
            var mesh = new THREE.Mesh(geometry, material);
        };
        // add mesh to planet meshes
        planetMeshes[planetName] = mesh;
        // Create and set orbit
        const orbit = createOrbit(orbitParams, planetData.renderParams.color, ORBIT_MESH_POINTS);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, 0, orbitParams.transformMatrix);
        mesh.position.set(pos.x, pos.y, pos.z);
        // Add to scene
        scene.add(orbit);
        scene.add(mesh);
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

    // Rotate Saturn's rings:
    // rings.rotation.z += 0.002; //rotate rings slightly

    controls.update();
    renderer.render(scene, camera);
}

// Start animation loop
requestAnimationFrame(animate);