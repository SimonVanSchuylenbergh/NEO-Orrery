// Imports
import * as THREE from 'https://cdn.skypack.dev/three@0.124.0/build/three.module.js';
// import { GLTFLoader } from 'https://cdn.jsdelivr.net/npm/three@0.114/examples/jsm/loaders/GLTFLoader.js';
import { OrbitControls } from 'https://cdn.skypack.dev/three@0.128.0/examples/jsm/controls/OrbitControls.js';
import { createOrbit, getOrbitPosition, JulianDateToTrueAnomaly } from './orbits.js'
import { JDToMJD, MJDToJD} from './TimeUtils.js'

// Constants
const DEG_TO_RAD = Math.PI / 180;
const TA_TIME_SCALE_FACTOR = 0.0001; // This will not be needed when the true anomaly code is included

const DEFAULT_MESH_N = 32;
const ORBIT_MESH_POINTS = 128;

const NEO_ORBIT_COLOR = 0x1e90FF;
const SHOWER_ORBIT_COLOR = 0x8B0000;
const PARENT_ORBIT_COLOR = 0xFF0000;

const NEO_COLOR = 0xFFFFFF;
const NEO_RADIUS = 0.01;
const MAX_VISIBLE_NEOS = 1;

const MOUSE_MIN_MOVE_CLICK = 0.005;

//starting time
let JD = (Date.now() / 86400000) + 2440587.5;
let MJD = JDToMJD(JD);

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

// Add lighting
const light = new THREE.AmbientLight(0x404040, 0.5); // Soft white light
scene.add(light);

// Setup Controls
let controls = null;

function resetCamera() {
    // resetting controls to remove damping --> https://codepen.io/boytchev/pen/qBJxdxo
    if (controls !== null) {controls.dispose( );}
    controls = new OrbitControls( camera, renderer.domElement );

    //camera.position.set(2, 2, 2);
    //camera.lookAt(0, 0, 0);
    controls.object.position.set(2, 2, 2);
    controls.target = new THREE.Vector3(0, 0, 0);
    controls.enableDamping = true;
}

resetCamera(); // Set camera position
controls.update(); //Is this needed?

// Setup Event Listeners
window.addEventListener("resize", () => {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);

    /*
    TESTING ONLY
    const canvas2d = document.getElementById('2DCanvas');
    const ctx = canvas2d.getContext('2d');
    ctx.canvas.width  = window.innerWidth;
    ctx.canvas.height = window.innerHeight;
    */
});

//All variables for key toggling go here

// all keyboard events are handled here
window.addEventListener("keydown", (event) => {

    if (event.code === 'Digit1') { //press to show a certain subset of the solar sytem bodies (e.g only planets or only high risk NEOs)

    }
    else if (event.code === 'Digit2') { //press to show a certain subset of the solar sytem bodies
        console.log('2 is pressed down!');
        
    }
    else if (event.code === 'Digit3') { //press to show a certain subset of the solar sytem bodies
        console.log('3 is pressed down!');
        
    }
    else if (event.code === 'Digit4') { //press to show a certain subset of the solar sytem bodies
        console.log('4 is pressed down!');
        
    }
    else if (event.code === 'KeyQ') { //decrease clicking threshold (logarithmically)
        console.log('Q key is pressed down!');
        
    }
    else if (event.code === 'KeyE') { //increase clicking threshold (logarithmically)
        console.log('E key is pressed down!');
        
    }
    else if (event.code === 'KeyA') { //cycle parameter in descending order (e.g. eccentricity)
        console.log('A key is pressed down!');
        
    }
    else if (event.code === 'KeyD') { //cycle parameter in ascending order
        console.log('D key is pressed down!');
        
    }
    else if (event.code === 'KeyW') { //reserved for whatever
        console.log('W key is pressed down!'); 
        
    }
    else if (event.code === 'KeyS') { //reserved for whatever
        console.log('S key is pressed down!');
    }
    else if (event.code === 'KeyT') { //track selected orbit
        if (highlightedObj !== null){
            controls.target = highlightedObj.userData.parent.bodyMesh.position;
        }
    }
    else if (event.code === 'KeyU') { //untrack selected orbit
        if (highlightedObj !== null){
            controls.target = highlightedObj.userData.parent.bodyMesh.position.clone();
        }
    }
    else if (event.code === 'Space') { //recenter on Sun
        controls.target = new THREE.Vector3(0, 0, 0);
    }
    else if (event.code === 'Backspace') { //recenter camera to initial settings
        resetCamera();
    }

});

const mouseDownXY = new THREE.Vector2(-10, -10);
const mouseUpXY = new THREE.Vector2(-10, -10);
const raycaster = new THREE.Raycaster(); //ray through the screen at the location of the mouse pointer (when the mouse is released)
raycaster.params.Line = {threshold: 0.05}; //this will just be a user interactive slider
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
        const allselectedObjects = raycaster.intersectObjects(scene.children);

        //remove duplicate intersections of the same orbit
        const seen = new Set();
        const selectedOrbits = allselectedObjects.filter(item => {
            if (!seen.has(item.object.uuid) && item.object.type === 'Line') {
                seen.add(item.object.uuid);
                return true; // unique uuid
            }
            return false; // same uuid
        });

        if (highlightedObj != null){ //clicking on the background deselects the current object (if there is one)
            highlightedObj.material.color.set(prevColor);
            highlightedObj = null;
            document.querySelector('.info-panel').style.display = "none";
        }
        
        if (selectedOrbits.length != 0){
            if (!moved) { stackedObjIndex = (stackedObjIndex + 1) % selectedOrbits.length; }
            
            highlightedObj = selectedOrbits[stackedObjIndex].object; //save the highlighted object
            prevColor = highlightedObj.material.color.getHex(); //save the highlighted object's previous color
            highlightedObj.material.color.set(0x00ff00); //highlight the select object if it is an orbit
            document.querySelector('.info-panel').style.display = "block";
            // Update info in the info panel
            const parentObj = selectedOrbits[stackedObjIndex].object.userData.parent;
            const obj_data = parentObj.data;

            document.getElementById('info-name').textContent = selectedOrbits[stackedObjIndex].object.userData.parent.name;
            document.getElementById('info-diameter').textContent = `Diameter: ${obj_data.extraParams.diameter} m`;
            document.getElementById('info-first-impact').textContent = `First possible impact: ${obj_data.extraParams.impact}`;
            document.getElementById('info-impact-period').textContent = `Possible impacts between ${obj_data.extraParams.years.split('-')[0]} and ${obj_data.extraParams.years.split('-')[1]}`;
            document.getElementById('info-risk').textContent = `Risk: ${obj_data.extraParams['PS max']}`;
            document.getElementById('info-vel').textContent = `Velocity: ${obj_data.extraParams.vel} km/s`;
            document.getElementById('info-a').textContent = `Semi-major axis: ${obj_data.orbitParams.a.toFixed(3)} AU`;
            document.getElementById('info-e').textContent = `Eccentricity: ${obj_data.orbitParams.e.toFixed(3)}`;
            document.getElementById('info-inc').textContent = `Inclination: ${(obj_data.orbitParams.inc / Math.PI * 180).toFixed(3)}\u00B0`;
            document.getElementById('info-node').textContent = `Longitude of ascending node: ${(obj_data.orbitParams.node / Math.PI * 180).toFixed(3)}\u00B0`;
            document.getElementById('info-peri').textContent = `Argument of perihelion: ${(obj_data.orbitParams.peri / Math.PI * 180).toFixed(3)}\u00B0`;
            document.getElementById('info-ma').textContent = `Mean anomaly: ${(obj_data.orbitParams.ma / Math.PI * 180).toFixed(3)}\u00B0`;
            document.getElementById('info-epoch').textContent = `Epoch: ${obj_data.orbitParams.epoch} (MJD)`;
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
    const sunTexture = sunTextureLoader.load('assets/body_textures/8k_sun.jpg');

    const geometry = new THREE.SphereGeometry(0.02, DEFAULT_MESH_N, DEFAULT_MESH_N);
    const material = new THREE.MeshBasicMaterial({map: sunTexture});
    const sunMesh = new THREE.Mesh(geometry, material);
    scene.add(sunMesh);
}

async function initializePlanets() {
    const planets_json = await readJSON('data/planet_data.json');
    for (const [planetName, planetData] of Object.entries(planets_json)) {
        const orbitParams = planetData.orbitParams;
        orbitParams.inc *= DEG_TO_RAD;
        orbitParams.node *= DEG_TO_RAD;
        orbitParams.peri *= DEG_TO_RAD;
        orbitParams.ma *= DEG_TO_RAD;
        // get planet texture
        const planetTextureName = planetData.renderParams.texture;
        const planetTextureLoader = new THREE.TextureLoader();
        const planetTexture = planetTextureLoader.load('assets/body_textures/' + planetTextureName);
        // check if planet is saturn's rings
        // if so, make it a ring geometry with specified parameters -- otherwise, make it a sphere egeometry
        // console.log(planetName)
        if (planetName == 'rings'){
            const geometry = new THREE.RingGeometry(planetData.renderParams.innerRadius, planetData.renderParams.outerRadius, 64);
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
        // Create and set orbit
        const orbit = createOrbit(orbitParams, planetData.renderParams.color, ORBIT_MESH_POINTS);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, 0, orbitParams.transformMatrix);
        mesh.position.set(pos.x, pos.y, pos.z);

        const body = new Body(planetName, planetData, orbit, mesh);

        orbit.userData.parent = body;
        mesh.userData.parent = body;
        planets.push(body);

        // Add to scene
        scene.add(orbit);
        scene.add(mesh);
    }
}

async function initializeNeos() {
    const neos_json = await readJSON('data/risk_list_neo_data.json');
    let i = 0;
    for (const [neoName, neoData] of Object.entries(neos_json)) {
        const orbitParams = neoData.orbitParams;
        orbitParams.inc *= DEG_TO_RAD;
        orbitParams.node *= DEG_TO_RAD;
        orbitParams.peri *= DEG_TO_RAD;
        orbitParams.ma *= DEG_TO_RAD;

        const geometry = new THREE.SphereGeometry(NEO_RADIUS, DEFAULT_MESH_N / 2, DEFAULT_MESH_N / 2);
        const material = new THREE.MeshBasicMaterial({ color: NEO_COLOR });
        const neoMesh = new THREE.Mesh(geometry, material);

        const orbit = createOrbit(orbitParams, NEO_ORBIT_COLOR, ORBIT_MESH_POINTS);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, 0, orbitParams.transformMatrix);
        neoMesh.position.set(pos.x, pos.y, pos.z);

        const body = new Body(neoName, neoData, orbit, neoMesh);

        orbit.userData.parent = body;
        neoMesh.userData.parent = body;
        planets.push(body);

        scene.add(orbit);
        scene.add(neoMesh);

        i += 1;
        if (i == MAX_VISIBLE_NEOS) { break };
    }
}

const animatedParentBodies = {}; // Store parent body positions and orbits

async function initializeShower() {
    let j = 0;
    let i = 0;
    const showers_json = await readJSON('data/stream_dataIAU2022.json');
    const parentBodies = await readJSON('data/stream_parentbody.json');
    const showerEntries = Object.entries(showers_json);
    for (const [showerName, showerData] of showerEntries) {
        const orbitParams = { ...showerData.orbitParams };
        const extraParams = showerData.extraParams;

        orbitParams.inc *= DEG_TO_RAD;
        orbitParams.node *= DEG_TO_RAD;
        orbitParams.peri *= DEG_TO_RAD;

        const orbit = createOrbit(orbitParams, SHOWER_ORBIT_COLOR, ORBIT_MESH_POINTS);
        scene.add(orbit);
        if (i + 1 < showerEntries.length) {
            const nextShowerData = showerEntries[i + 1][1];
            if (showerData.extraParams.Code !== nextShowerData.extraParams.Code) {
                j += 1;
                for (const [parentBodyName, parentBodyData] of Object.entries(parentBodies)) {
                    if (parentBodyData.extraParams.Code === extraParams.Code) {
                        const orbitParams_parent = parentBodyData.orbitParams;
                        orbitParams_parent.inc *= DEG_TO_RAD;
                        orbitParams_parent.node *= DEG_TO_RAD;
                        orbitParams_parent.peri *= DEG_TO_RAD;
                        orbitParams_parent.ma *= DEG_TO_RAD;

                        const geometry = new THREE.SphereGeometry(NEO_RADIUS, DEFAULT_MESH_N / 2, DEFAULT_MESH_N / 2);
                        const material = new THREE.MeshBasicMaterial({ color: PARENT_ORBIT_COLOR });
                        const parentMesh = new THREE.Mesh(geometry, material);

                        const orbit = createOrbit(orbitParams_parent, PARENT_ORBIT_COLOR, ORBIT_MESH_POINTS);
                        const pos = getOrbitPosition(orbitParams_parent.a, orbitParams_parent.e, 0, orbitParams_parent.transformMatrix);
                        parentMesh.position.set(pos.x, pos.y, pos.z);

                        const body = new Body(parentBodyName, parentBodyData, orbit, parentMesh)
                        parentMesh.userData.parent = body;
                        orbit.userData.parent = body;
                        neos.push(body)

                        scene.add(orbit);
                        scene.add(parentMesh);

                        animatedParentBodies[parentBodyName] = {
                            mesh: parentMesh,
                            orbitParams: orbitParams_parent,
                            currentAnomaly: orbitParams_parent.ma
                        };

                        break;
                    }
                }
            }
        }
        i += 1;
        if (j === MAX_VISIBLE_NEOS) {
            break;
        }
    }
}

function updateParentBodyPosition(parentBody, JD) {
    const orbitParams = parentBody.orbitParams;
    const trueAnomaly = JulianDateToTrueAnomaly(orbitParams, JD);

    const pos = getOrbitPosition(orbitParams.a, orbitParams.e, trueAnomaly, orbitParams.transformMatrix);
    parentBody.mesh.position.set(pos.x, pos.y, pos.z);
}

// Add radial gradient plane
function createRadialGradientPlane(width, height) {
    const geometry = new THREE.PlaneGeometry(width, height, 1, 1);
    const material = new THREE.ShaderMaterial({
        vertexShader: `
            varying vec2 vUv;
            void main() {
                vUv = uv;
                gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }
        `,
        fragmentShader: `
            varying vec2 vUv;
            void main() {
                float distanceFromCenter = length(vUv - vec2(0.5, 0.5));
                float alpha = (1.0 - distanceFromCenter * 2.0)*0.5;
                alpha = clamp(alpha, 0.0, 1.0);
                if (alpha < 0.01) {
                    discard;
                }
                gl_FragColor = vec4(1.0, 0.0, 0.0, alpha);
            }
        `,
        transparent: true,
        side: THREE.DoubleSide,
        depthWrite: false,
        depthTest: false,
    });

    const plane = new THREE.Mesh(geometry, material);
    plane.rotation.x = Math.PI / 2;
    plane.renderOrder = 0;
    return plane;
}

// Function to create a radial gradient plane with an exponential drop-off
function createSunGradientPlane(width, height) {
    const geometry = new THREE.PlaneGeometry(width, height, 1, 1);
    const material = new THREE.ShaderMaterial({
        vertexShader: `
            varying vec2 vUv;
            void main() {
                vUv = uv;
                gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
            }
        `,
        fragmentShader: `
            varying vec2 vUv;
            void main() {
                // Calculate distance from the center of the plane (0.5, 0.5) in UV space
                float distanceFromCenter = length(vUv - vec2(0.5, 0.5));

                // Exponential drop-off for the intensity
                float alpha = exp(-50.0 * distanceFromCenter);  // Adjust for a faster/slower fade

                // Clamp alpha to ensure it's between 0 and 1
                alpha = clamp(alpha, 0.0, 1.0);

                // Discard very transparent fragments
                if (alpha < 0.01) {
                    discard;
                }

                // Set the color to a white-yellowish tone with the calculated alpha
                gl_FragColor = vec4(1.0, 0.95, 0.6, alpha);
            }
        `,
        transparent: true,
        side: THREE.DoubleSide,
        depthWrite: false,
        depthTest: false,
    });

    const plane = new THREE.Mesh(geometry, material);
    plane.renderOrder = 1;  // Ensure the plane renders after other objects
    return plane;
}

// Function to create a single plane that always faces the camera
function createBillboardPlane(size) {
    const width = size;
    const height = size;

    // Create the radial gradient plane
    const gradientPlane = createSunGradientPlane(width, height);

    // Return the gradient plane mesh
    return gradientPlane;
}

// Example usage: Create a billboard plane that follows the camera
const hazeSize = 2.0;  // Set the size of the haze plane
const billboardPlane = createBillboardPlane(hazeSize);
scene.add(billboardPlane);

// Function to update the plane to always face the camera (billboarding effect)
function updateBillboard(plane, camera) {
    // Set the plane's rotation to always face the camera
    plane.lookAt(camera.position);
}


const planeWidth = 5.204 * 2;
const radialGradientPlane = createRadialGradientPlane(planeWidth, planeWidth);
// scene.add(radialGradientPlane);


// Data
let sunMesh;
const planets = [];
const neos = [];


class Body {
    constructor(name, data, orbitMesh, bodyMesh) {
        this.name = name;
        this.data = data;
        this.orbitMesh = orbitMesh;
        this.bodyMesh = bodyMesh;

        // These will be pointers used for cycling through objects
        this.nextLargerSize;
        this.nextSmallerSize;
        this.nextLargerRisk;
        this.nextSmallerRisk;
        this.nextLargerImpactTime;
        this.nextSmallerImpactTime;
        this.nextLargerA;
        this.nextSmallerA;
        this.nextLargerE;
        this.nextSmallerE;
    }

    setPosition(pos) {
        this.bodyMesh.position.set(pos.x, pos.y, pos.z)
    }
}



addSun();
await initializePlanets(); // Initialize planets once
await initializeNeos(); // Initialize NEOs once
await initializeShower();
// console.log(planets.Saturn);


/*
//FOR TESTING PURPOSES
// Function to draw nested ellipses
function drawNestedEllipses(x, y, initialWidth, initialHeight, count) {
    const canvas = document.getElementById('2DCanvas');
    const ctx = canvas.getContext('2d');
    ctx.canvas.width  = window.innerWidth;
    ctx.canvas.height = window.innerHeight;


    for (let i = 0; i < count; i++) {
        ctx.beginPath();
        ctx.ellipse(
            x,                  // Center x
            y,                  // Center y
            initialWidth - i * 10,  // Semi-major axis (width)
            initialHeight - i * 10,  // Semi-minor axis (height)
            0.3,                  // Rotation
            0, Math.PI * 2     // Start and end angle
        );
        ctx.strokeStyle = "white";
        ctx.stroke();
    }
}

drawNestedEllipses(window.innerWidth / 2, window.innerHeight / 2, 150, 100, 10); // x, y, initial width, initial height, number of ellipses
*/

// Animation loop with FPS control
function animate(time) {
    requestAnimationFrame(animate);

    // Limit frame rate
    const deltaTime = time - lastFrameTime;
    if (deltaTime < frameInterval) {
        return; // Skip frame if too soon
    }
    lastFrameTime = time;

    let timeSpeed = 864000;
    let deltaJulian = deltaTime * timeSpeed / 1000 / 86400;

    JD += deltaJulian;
    MJD += deltaJulian;

    // Update planet positions and rotation
    for (let i = 0; i < planets.length; i++) {
        const orbitParams = planets[i].data.orbitParams;
        const extraParams = planets[i].data.extraParams;
        const trueAnomaly = JulianDateToTrueAnomaly(orbitParams, JD);
        // Update Position
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, trueAnomaly, orbitParams.transformMatrix);
        planets[i].setPosition(pos);
        // Rotate
        //planets[i].bodyMesh.rotation.x += orbitParams.rotateX;
        //planets[i].bodyMesh.rotation.y += orbitParams.rotateY;
        //planets[i].bodyMesh.rotation.z += orbitParams.rotateZ;
    }

    // Update NEO positions
    for (let i = 0; i < neos.length; i++) {
        const orbitParams = neos[i].data.orbitParams;
        // console.log(neos[i])
        const trueAnomaly = JulianDateToTrueAnomaly(orbitParams, MJD);
        const pos = getOrbitPosition(orbitParams.a, orbitParams.e, trueAnomaly, orbitParams.transformMatrix);
        neos[i].setPosition(pos);
    }

    for (const parentBodyName in animatedParentBodies) {
        updateParentBodyPosition(animatedParentBodies[parentBodyName], JD);
    }

    // Update the billboard plane to face the camera
    updateBillboard(billboardPlane, camera);

    controls.update();
    renderer.render(scene, camera);
}

// Start animation loop
requestAnimationFrame(animate);