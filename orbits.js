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