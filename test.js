// let planets = fetch('data/planet_data.json'); // import planets from json
// let neos = fetch('data/risk_list_neo_data.json'); // import neos from json

let planets = require('./data/planet_data.json') 
let neos = require('./data/risk_list_neo_data.json')

console.log(planets.Mercury)
console.log(neos['2005QK76'])