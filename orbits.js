import * as THREE from 'https://cdn.skypack.dev/three@0.128.0/build/three.module.js';

export function createOrbit(orbitParams, color, n_mesh_points) {
    const orbit_segment_const = 2 * Math.PI / n_mesh_points;
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

    for (let i = 0; i <= n_mesh_points; i++) {
        const eccentric_anomaly = orbit_segment_const * i; // Angle
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

export function getOrbitPosition(a, e, trueAnomaly, matrix) {
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

export function JulianDateToTrueAnomaly(orbitParams, JD) {
    const newMA = getCurrentMeanAnomaly(orbitParams.a, orbitParams.ma, JD, orbitParams.epoch);
    const E = solveKepler(orbitParams.e, newMA);
    return computeTrueAnomaly(E, orbitParams.e);
}

function getCurrentMeanAnomaly(a, ma, JD, epoch) { //JD and epoch need to use the same reference system (either MJD or JD) 
    const mu = 0.0002959122082855911025;
    return (JD - epoch) * Math.sqrt(mu / Math.abs(a**3)) + ma;
}

//Computes the the true anomaly, given the eccentric anomaly along with the eccentricity
function computeTrueAnomaly(E, e) { return 2*Math.atan(Math.sqrt((1+e) / (1-e)) * Math.tan(E/2)) }

//function for solving Kepler's equation using a binary search approach given the eccentricity and the mean anomaly
function solveKepler(e, M) {
    const espLim = 10*Math.max(Number.EPSILON, Math.abs(M)*Number.EPSILON);

    if (e == 0) { return M } //trivial case

    const keplerFunc = (e, E) => { return E - e * Math.sin(E); };
    let E = M; // starting guess
    const EMult = Math.sqrt(2);
        
    let minBound = 0;
    let maxBound = 0;

    let MTest = 0;
    let MDiff = 0;

    //initialize min and max bounds
    if (M < 0) {
        while (true) {
            MTest = keplerFunc(e, E);
            MDiff = M - MTest;
            if (Math.abs(MDiff) < espLim) { return E; }
            if (MDiff > 0) {
                minBound = E;
                break;
            }
            else { 
                maxBound = E;
                E *= EMult;
            }
        }
    }
    else {
        while (true) {
            MTest = keplerFunc(e, E);
            MDiff = M - MTest;
            if (Math.abs(MDiff) < espLim) { return E; }
            if (MDiff > 0) {
                minBound = E;
                E *= EMult;
            }
            else { 
                maxBound = E;
                break;
            }
        }
    }

    //let i = 0;
    //let maxI = 20;

    while (true) { //perform a binary search to solve for E
        E = (maxBound + minBound) / 2; //take E at the midpoint
        MTest = keplerFunc(e, E);
        MDiff = M - MTest;
        if (Math.abs(MDiff) < espLim) { return E; }
        if (MDiff > 0) { minBound = E; }
        else { maxBound = E; }
        //i++;
        //console.log(E, MTest, M, MDiff);
        //if (i == maxI) { break; }
    }
}
