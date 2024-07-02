import { ConstraintDistance, ConstraintPin } from "./constraints.js";
import { translate, rotate, getNeighborsDelauney, calculateInertiaMoment,
    integrateOde, integrateOdeOnTarray, isShapeEqual} from "./math_utils.js";
//IntegrateSystem(fun, x0, tf, h) where x0 is a R dim array, fun is a R->R dxdt fun 
// integrateOdeOnTarray(fun,xo,tArray,method = "RK4")
//h is the time step, tf is the final time, by default uses RK4
//all constraints must have getConstraintValue(xarray) getJacobianParticles(xarray) getDotJacobianParticles(varray) methods
//function dArraydtFun(array, h){
    //calculates numericalderivative for an array of values, assuming they come as discrete evaluation of a t->R function sampled at fixed intervals h
    //returns an array of the same length as array
const DEBUG = true;

if (DEBUG){
    console.log("solver.js loaded","PENDING TESTING SPRINGS AND DAMPING FORCES");
}






/*
GENERAL FUNCTIONALITIES
*/
//basic function to calculate COM given positions array xs [Nparticles, 2] and masses

function calculateCOM(xs,masses){
    //xs array of particle possitions
     if (xs.length !== masses.length){
        throw new Error("The xs and masses arrays should have the same length")
      }
    if (xs[0].length !== 2){
      throw new Error("The xs array should have 2 components")
    }
    let total_mass = masses.reduce((a,b)=>a+b,0);
    let weighted_xs = xs.map((x,i)=>math.multiply(masses[i]/total_mass,x));
    let COM = weighted_xs.reduce((a,b)=>math.add(a,b),math.zeros(2));
  
      
      return COM;
    }



function calculateKineticEnergy(vs,masses){
    if (vs.length !== masses.length){
        throw new Error("The xs and masses arrays should have the same length")
    }
    if (vs[0].length !== 2){
    throw new Error("The xs array should have 2 components")
    }
    //Ekin = 0.5*v^{T}Mv
    let M = math.matrix(masses.map(mass=>[mass,mass])).reshape([-1]);//vx and vy have the same mass(same particle)
    M = math.diag(M);
    let v_matrix = math.flatten(math.matrix(vs));

    let Ekin = 0.5*math.multiply(math.transpose(v_matrix),math.multiply(M,v_matrix));
    return Ekin;
    }










/*
NON-CONTACT FORCES CLASSES, signature is simple, in the constructor we pass the particle indices and extra specific force parameters
All classes should have a getForceArray method that returns a [nparticles,2] array with the forces acting on each particle
The input to getForceArray should be the current positions, velocities and masses of the particles
*/

class ConstantForce{
    constructor(particle_indices, force){
        this.particle_indices = particle_indices;
        this.force = force;
    }
    getForceArray(xarray, varray, masses){
        if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
        throw new Error("xarray, varray and masses should be arrays")
        }
        let nparticles = xarray.length;
        let nforces = this.particle_indices.length;
        let force_matrix= math.zeros(nparticles,2,"sparse");
        for (let i = 0; i < nforces; i++){
        let particle_index = this.particle_indices[i];
        let force = this.force;
        let force_x = force[0];
        let force_y = force[1];
        force_matrix.set([particle_index, 0], force_x);
        force_matrix.set([particle_index, 1], force_y);
        }
        return force_matrix.toArray();
    }
    }

class Gravity{
    constructor(particle_indices,g = 9.8,direction = [0,-1]){
        this.g = g;
        this.particle_indices = particle_indices;
        this.direction = direction;
    }
    getForceArray(xarray, varray, masses){
        //lets throw error if xarray, varray and masses are not arrays
        if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){

        throw new Error("xarray, varray and masses should be arrays")
        }
        let nparticles = xarray.length;
        let nforces = this.particle_indices.length;
        let force_matrix = math.zeros(nparticles,2,"sparse");
        for (let i = 0; i < nforces; i++){
        let particle_index = this.particle_indices[i];
        let mass = masses[particle_index];
        let force = this.g*mass;
        let force_x = force * this.direction[0];
        let force_y = force * this.direction[1];
        force_matrix.set([particle_index, 0], force_x);
        force_matrix.set([particle_index, 1], force_y);
        }
        return force_matrix.toArray();
    }
}


class Damping{
    constructor(particle_indices, damping_coefficient){
        this.particle_indices = particle_indices;
        this.damping_coefficient = damping_coefficient;
    }
    getForceArray(xarray, varray, masses){
        if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
        throw new Error("xarray, varray and masses should be arrays")
        }
            let nparticles = xarray.length;
            let nforces = this.particle_indices.length;
            let force_matrix = math.zeros(nparticles,2,"sparse");
            for (let i = 0; i < nforces; i++){
            let particle_index = this.particle_indices[i];
            let damping_coefficient = this.damping_coefficient;
            let velocity = varray[particle_index];
            let force_x = -damping_coefficient*velocity[0];
            let force_y = -damping_coefficient*velocity[1];
            force_matrix.set([particle_index, 0], force_x);
            force_matrix.set([particle_index, 1], force_y);
        }
        return force_matrix.toArray();
    }
}

class Spring{
    constructor(idx1,idx2, spring_constant, rest_length){
        this.particle_indices = [idx1,idx2];
        this.spring_constant = spring_constant;
        this.rest_length = rest_length;
    }
    getForceArray(xarray, varray, masses){
        if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
        throw new Error("xarray, varray and masses should be arrays")
        }
        let nparticles = xarray.length;
        let nforces = this.particle_indices.length;
        let force_matrix = math.zeros(nparticles,2,"sparse");
        let idx1 = this.particle_indices[0];
        let idx2 = this.particle_indices[1];
        let x1 = xarray[idx1];
        let x2 = xarray[idx2];
        let x1_x = x1[0];
        let x1_y = x1[1];
        let x2_x = x2[0];
        let x2_y = x2[1];
        let distance = Math.sqrt((x2_x-x1_x)**2 + (x2_y-x1_y)**2);
        let force_magnitude = this.spring_constant*(distance-this.rest_length);
        let force_x = force_magnitude*(x2_x-x1_x)/distance;
        let force_y = force_magnitude*(x2_y-x1_y)/distance;
        force_matrix.set([idx1, 0], force_x);
        force_matrix.set([idx1, 1], force_y);
        force_matrix.set([idx2, 0], -force_x);
        force_matrix.set([idx2, 1], -force_y);
        return force_matrix.toArray();
    }
}

class Spring2Point{
    //This spring acts between a target particle and a specified fixed point,this fixed point can be set with mouse coordinates later
    constructor(idx1, fixed_point, spring_constant, rest_length){
        this.particle_indices = [idx1];
        this.spring_constant = spring_constant;
        this.rest_length = rest_length;
        this.fixed_point = fixed_point;
    }
    getForceArray(xarray, varray, masses){
        if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
        throw new Error("xarray, varray and masses should be arrays")
        }
        let nparticles = xarray.length;
        let nforces = this.particle_indices.length;
        let force_matrix = math.zeros(nparticles,2,"sparse");
        let idx1 = this.particle_indices[0];
        let x1 = xarray[idx1];
        let x1_x = x1[0];
        let x1_y = x1[1];
        let x2_x = this.fixed_point[0];
        let x2_y = this.fixed_point[1];
        let distance = Math.sqrt((x2_x-x1_x)**2 + (x2_y-x1_y)**2);
        let force_magnitude = this.spring_constant*(distance-this.rest_length);
        let force_x = force_magnitude*(x2_x-x1_x)/distance;
        let force_y = force_magnitude*(x2_y-x1_y)/distance;
        force_matrix.set([idx1, 0], force_x);
        force_matrix.set([idx1, 1], force_y);
        return force_matrix.toArray();
    }
}



function computeExternalForces(forces,xarray,varray,masses){
    /*
    * forces is an array of force objects
    * This returns a [nparticles,2] array with the external forces acting on each particle
    */
    let nparticles = xarray.length;
    let external_forces = math.zeros(nparticles,2,"sparse");
    for (let force of forces){
      let ext_force = force.getForceArray(xarray, varray, masses);
  
      external_forces = math.add(external_forces, ext_force); 
    }
    return external_forces;
  }










/*
SOLVE CONSTRAINTS FUNCTIONS (NON-CONTACT)
*/
function computeJacobiansNonContact(constraints,xarray,varray){
    /*
    *   Compute the jacobian matrix J and dot jacobian matrix dot_J for the constraints
    *   J is a matrix of size [nconstraints, 2*nparticles]
    *  dot_J is a matrix of size [nconstraints, 2*nparticles]
    *  The constraints should be an array of constraint objects
    */

    // type checking, constraints,xarray,varray should be arrays
    if (!Array.isArray(constraints) || !Array.isArray(xarray) || !Array.isArray(varray)){
        throw new Error("constraints, xarray and varray should be arrays")
    }
    const nconstraints = constraints.length;
    const nparticles = xarray.length;
    const J = math.zeros(nconstraints, 2*nparticles,"sparse");
    const dot_J = math.zeros(nconstraints, 2*nparticles,"sparse");

    for (let i = 0; i < nconstraints; i++){
        let constraint = constraints[i];
        let constraint_particle_indices = constraint.getParticleIndices();
        let jacobian_particles = constraint.getJacobianParticles(xarray);
        let dot_jacobian_particles = constraint.getDotJacobianParticles(varray);

        if (constraint_particle_indices.length !== jacobian_particles.length){
        throw new Error("The number of particles in the constraint is not equal to the number of jacobian particles")
        }
        if (constraint_particle_indices.length !== dot_jacobian_particles.length){
        throw new Error("The number of particles in the constraint is not equal to the number of dot jacobian particles")
        }


        for (let j = 0; j < constraint_particle_indices.length; j++){
        let particle_index = constraint_particle_indices[j];
        let jacobian_particle = jacobian_particles[j];
        if (jacobian_particle.length !== 2){
            throw new Error("The jacobian particle should have 2 components")
        }
        let jacobian_particle_x = jacobian_particle[0];
        let jacobian_particle_y = jacobian_particle[1];
        J.set([i, 2*particle_index], jacobian_particle_x);
        J.set([i, 2*particle_index + 1], jacobian_particle_y);

        let dot_jacobian_particle = dot_jacobian_particles[j];
        if (dot_jacobian_particle.length !== 2){
            throw new Error("The dot jacobian particle should have 2 components")
        }
        let dot_jacobian_particle_x = dot_jacobian_particle[0];
        let dot_jacobian_particle_y = dot_jacobian_particle[1];
        dot_J.set([i, 2*particle_index], dot_jacobian_particle_x);
        dot_J.set([i, 2*particle_index + 1], dot_jacobian_particle_y);

        
        }
    }

    return [J, dot_J];
    }


 //lets create the non dummy solve_constraints function
  //With Baumgarte stabilization the differentiated constraint eq is dotJv+Ja+alpha*dotC+beta*C = 0
  // Replacing a = M^{-1}(Fext + Fconstraint) and dotC = Jv we get dotJv + J(M^{-1}(Fext + Fconstraint)) + alpha*Jv + beta*C = 0
//And Fconstraint = J^T*lambda, so we get dotJv + J(M^{-1}(Fext + J^T*lambda)) + alpha*Jv + beta*C = 0
// so the resulting linear system to solve in A*lambda = b  where A = J*M^{-1}*J^T  and 
//b = -dotJv - J*M^{-1}*Fext - alpha*Jv - beta*C
// lets implement this
function solveNonContactConstraints(J, dot_J, xarray, varray,
    external_forces, masses, constraints,
     alpha = 0.00, beta = 0.00){
    /*
    *  ONLY FOR NON-CONTACT CONSTRAINTS
    */

    //lets do type checking, J,dot_J external forces mathjs matrices, xarray, varray, masses arrays
    if (!math.isMatrix(J) || !math.isMatrix(dot_J) || !math.isMatrix(external_forces)){
    throw new Error("J, dot_J and external_forces should be mathjs matrices")
    }
    if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
    throw new Error("xarray, varray and masses should be arrays")
    }
    let nparticles = xarray.length;
    let nconstraints = J.size()[0];
    let M_inv = math.matrix(masses.map(mass=>[1/mass,1/mass])).reshape([-1]);//vx and vy have the same mass(same particle)
    M_inv = math.diag(M_inv,"sparse");

    let A = math.multiply(J, math.multiply(M_inv, math.transpose(J)));
    //lets check shape of A, it should be [nconstraints,nconstraints]
    if (A.size()[0] !== nconstraints || A.size()[1] !== nconstraints){
    throw new Error("The A matrix should have the same number of constraints")
    }
    //lets compute b = -dot_J*varray - J*M_inv*external_forces - alpha*J*varray - beta*constraints
    //Lets do it in parts and then sum
    //first varray [nparticles,2] to matrix and flatten-> [2*nparticles]
    //the same for external forces [nparticles,2] to matrix and flatten-> [2*nparticles]
    //Get constraints_vals [nconstraints] from constraints

    let constraints_vals = math.matrix(
        constraints.map(constraint => constraint.getConstraintValue(xarray))
        );
    let vmatrix = math.flatten(varray);
    external_forces =math.flatten(external_forces);
    let part1 = math.multiply(dot_J, vmatrix);
    let part2 = math.multiply(J, math.multiply(M_inv, external_forces));
    let part3 = math.multiply(alpha, math.multiply(J, vmatrix));
    let part4 = math.multiply(beta, constraints_vals).reshape([nconstraints,1]);


    let b = math.add(math.add(part1, part2), math.add(part3, part4));
    //all parts appear with a minus sign in the equation
    b = math.multiply(b,-1);
    let lagrange_multipliers = math.lusolve(A,b).reshape([-1]);
    // lets create a dict-like object to store the results
    let out_dict = {
    lagrange_multipliers: lagrange_multipliers,
    constraints_vals: constraints_vals,
    }
    return out_dict;

}
















/*
TIME STEPPERS
*/

function stepTestSemiEuler(STATE, PARAMETERS){
    /*
    *   Provisional semi-implicit euler integration
    */
    let xs = STATE.xs;
    let vs = STATE.vs;
    let masses = STATE.masses;
    let external_forces = STATE.external_forces;
    let constraint_forces = STATE.constraint_forces;
    let contact_forces = STATE.contact_forces;
    let friction_forces = STATE.friction_forces;
    let total_forces = xs.map((x,i) => math.add(external_forces[i], constraint_forces[i], contact_forces[i], friction_forces[i]));

    let dt = PARAMETERS["dt"];
    // gravity is assumed to be computed already in the external forces

    let new_vs = vs.map((v,i) => math.add(v, math.multiply(total_forces[i], dt/masses[i])));
    let new_xs = xs.map((x,i) => math.add(x, math.multiply(new_vs[i], dt)));
    let new_time = STATE.time + dt;


    let newState = {
        "xs": new_xs,
        "vs": new_vs,
        "time": new_time,
    }

    // lets add the rest of properties to newState copying them from STATE, excluding ["xs", "vs","time"]
    for (let key in STATE){
        if (!["xs", "vs", "time"].includes(key)){
            newState[key] = STATE[key];
        }
    }

    return newState;
}


// lets write a semi-implicit euler time stepper that takes into account constraints
// we'll skip so by now the contact forces and friction forces resolution


function initState(STATE, PARAMETERS){
    
    //we need to create external_forces, constraint_forces, contact_forces, friction_forces , we'll initialize them to zeros    

    if (!PARAMETERS.hasOwnProperty("alpha") || !PARAMETERS.hasOwnProperty("beta")){
        throw new Error("PARAMETERS must have keys alpha and beta")
    }

    // lets throw an error if STATE has not the required minimum keys, xs, vs, masses
    let mandatory = ["xs", "vs", "masses"]
    for (let key of mandatory){
        if (!STATE.hasOwnProperty(key)){
            throw new Error(`The STATE object must have the key ${key}`)
        }
    }

    let [alpha, beta] = [PARAMETERS.alpha, PARAMETERS.beta];
    let s = STATE;
    let nparticles = s.xs.length;
    let [constraintsObjs, forceObjs] = stepUnpackState(s);
    let externalForcesArray = computeExternalForces(forceObjs, s.xs, s.vs, s.masses);
    let [J, dot_J] = computeJacobiansNonContact(constraintsObjs, s.xs, s.vs);
    let outSolveConstraints = solveNonContactConstraints(J, dot_J, s.xs, s.vs,
                                                            externalForcesArray, s.masses,
                                                             constraintsObjs, alpha, beta);
    let constraintForces= math.multiply(math.transpose(J), outSolveConstraints.lagrange_multipliers).reshape([nparticles,2]);
    let frictionForces = math.zeros(nparticles,2);
    let contactForces = math.zeros(nparticles,2);
    let totalForces = math.add(externalForcesArray, constraintForces, contactForces, frictionForces);

    s["external_forces"] = externalForcesArray.toArray();
    s["constraint_forces"] = constraintForces.toArray();
    s["contact_forces"] =   contactForces.toArray();
    s["friction_forces"] = frictionForces.toArray();
    s["total_forces"] = totalForces.toArray();
    s["constraints_vals"] = outSolveConstraints.constraints_vals.toArray();
    s["lagrange_multipliers"] = outSolveConstraints.lagrange_multipliers.toArray();
    s["J"] = J.toArray();
    s["dot_J"] = dot_J.toArray();
    s["time"] = 0;

    return s;
}


function stepUnpackState(STATE){
    // We turn the constraint info into constraint objects
    // additionally we turn force info into force objects
    // I'll enumerate the field and format  of the different entries in STATE, relevant for this function
    //constraints_distance: [[idx1, idx2, distance],...] ConstraintDistance class has signature  constructor(particle1_idx, particle2_idx, distance)
    // constraints_pin: [[idx,pin_point, pin_distance],...] ConstraintPin class has signature constructor(particle_idx, pin_point, pin_distance)
    // constant_forces: [[idx1, force1],...] 
    // gravity: [g, direction] Gravity class has signature constructor(particle_indices, g, direction)
    // by default particles_indices arreys are [0,1,2,...,Nparticles-1]
    // damping: [[idx1, damping_coefficient],...] Damping class has signature constructor(particle_indices, damping_coefficient)
    // springs: [[idx1,idx2, spring_constant, rest_length],...] Spring class has signature constructor(idx1,idx2, spring_constant, rest_length)

    //lets throw an error if the following keys are not present
    let mandatory = ["constraints_distance", "constraints_pin", "xs", "vs", "masses"]
    for (let key of mandatory){
        if (!STATE.hasOwnProperty(key)){
            throw new Error(`The STATE object must have the key ${key}`)
        }
    }

    let s = STATE;

    let constraintsDistanceObjs = s.constraints_distance.map(c => new ConstraintDistance(c[0], c[1], c[2]));
    let constraintsPinObjs = s.constraints_pin.map(c => new ConstraintPin(c[0], c[1], c[2]));
    let constraintsObjs = constraintsDistanceObjs.concat(constraintsPinObjs);
    let nparticles = s.xs.length;
    let nparticles_array = Array.from(Array(nparticles).keys());
    let forceObjs = [];
    //if constant_forces/gravity/damping/springs not present, we dont add nothing to forceObjs
    if (s.hasOwnProperty("constant_forces")){
        forceObjs = s.constant_forces.map(f => new ConstantForce([f[0]], f[1]));
    }
    if (s.hasOwnProperty("gravity")){
        forceObjs.push(new Gravity(nparticles_array, s.gravity[0], s.gravity[1]));
    }
    if (s.hasOwnProperty("damping")){
        forceObjs = forceObjs.concat(s.damping.map(f => new Damping([f[0]], f[1])));
    }
    if (s.hasOwnProperty("springs")){
        forceObjs = forceObjs.concat(s.springs.map(f => new Spring(f[0], f[1], f[2], f[3])));
    }



    return [constraintsObjs, forceObjs];

}

function stepSemiEulerNonContact(STATE, PARAMETERS){
    // Lets throw an error if STATE and PARAMETERS doesnt have the required keys
    // STATE must have, well on second thought, lets pack the STATE preprocessing in a function

    //throw an error if alpha/beta not in PARAMETERS
    if (!PARAMETERS.hasOwnProperty("alpha") || !PARAMETERS.hasOwnProperty("beta")){
        throw new Error("PARAMETERS must have keys alpha and beta")
    }

    //throw an error if dt not in PARAMETERS
    if (!PARAMETERS.hasOwnProperty("dt")){
        throw new Error("PARAMETERS must have key dt")
    }
    let s = STATE;
    let nparticles = s.xs.length;
    let dt = PARAMETERS["dt"];
    let [alpha, beta] = [PARAMETERS.alpha, PARAMETERS.beta];
    let [constraintsObjs, forceObjs] = stepUnpackState(s);
    let externalForcesArray = computeExternalForces(forceObjs, s.xs, s.vs, s.masses);
    let [J, dot_J] = computeJacobiansNonContact(constraintsObjs, s.xs, s.vs);
    let outSolveConstraints = solveNonContactConstraints(J, dot_J, s.xs, s.vs,
                                                            externalForcesArray, s.masses,
                                                             constraintsObjs, alpha, beta);

    let lagrange_multipliers = outSolveConstraints.lagrange_multipliers;
    let constraintVals = outSolveConstraints.constraints_vals;
    let constraintForces= math.multiply(math.transpose(J), lagrange_multipliers).reshape([nparticles,2]);


    if (constraintForces.size()[0] !== nparticles){
        throw new Error("The constraint forces should have the same number of particles")
      }

    let totalForces= math.add(externalForcesArray, constraintForces);

    
    let accs = math.zeros(nparticles,2);// acceleration matrix
    for (let i = 0; i < nparticles; i++){
      let mass = s.masses[i];
      let forceX = totalForces.get([i,0]);
      let aX = forceX / mass;
      let forceY = totalForces.get([i,1]);
      let aY = forceY / mass;
      accs.set([i,0], aX);
      accs.set([i,1], aY);
    }

    let vMatrix = math.matrix(s.vs);
    let xMatrix = math.matrix(s.xs);
    let newVmatrix = math.add(vMatrix, math.multiply(accs, dt));
    let newXarray = math.add(xMatrix, math.multiply(newVmatrix, dt));
    let newTime = s.time + dt;

    //now lets create the new state object (everything that is matrix is converted to array)
    //The new state will have
    //xs,vs,accs,external_forces,constraint_forces,contact_forces and friction_forces (zeros(nparticles,2) arrays)
    //constraints_vals, lagrange_multipliers, time
    // J, Jdot 
    let newState = {
        "J": J.toArray(),
        "dot_J": dot_J.toArray(),
        "xs": newXarray.toArray(),
        "vs": newVmatrix.toArray(),
        "accs": accs.toArray(),
        "external_forces": externalForcesArray.toArray(),
        "constraint_forces": constraintForces.toArray(),
        "contact_forces": math.zeros(nparticles,2).toArray(),
        "friction_forces": math.zeros(nparticles,2).toArray(),
        "total_forces": totalForces.toArray(),//this is the sum of all forces acting on the particles
        "constraints_vals": constraintVals.toArray(),
        "lagrange_multipliers": lagrange_multipliers.toArray(),
        "time": newTime,
    }  

    //for the rest of the properties ( that are not in newState) we'll copy them from the original STATE object

    let keysNewState = Object.keys(newState);
    for (let key in s){
        if (!keysNewState.includes(key)){
            newState[key] = s[key];
        }
    }
    return newState;

}
   
    










/*
CONTACT AND COLLISIONS
*/


class Polygon{
    /* Basic usage of the Polygon class
    * let vertices = [[0,0],[1,0],[1,1],[0,1]];
    * let polygon = new Polygon(vertices, 0);
    * let point = [0.5,0.5];
    * let [closest_edge,distance] = polygon.getClosestEdge(point);
    * let closest_point = polygon.getClosestPoint2Edge(point,polygon.edges[closest_edge]);
    * let normal = polygon.getNormal(closest_edge);
    * let normal_projection = polygon.getClosestEdgeNormalProjection(point);
    * let normal_projection_vector = polygon.getClosestEdgeNormalProjectionVector(point);
    * let distance2edge = polygon.getDistance2Edge(point,polygon.edges[closest_edge]);
    * let center = polygon.getCenter();
    * let edge_midpoint = polygon.getEdgeMidpoint(closest_edge);
    * let normals = polygon.getNormals();
    * let closest_edge_normal_projection = polygon.getClosestEdgeNormalProjection(point);
    */
    constructor(vertices, id){
        this.vertices = vertices;
        this.edges = this.getEdges();
        this.id = id;
    }

    updateVertices(vertices){
        this.vertices = vertices;
        this.edges = this.getEdges();

    }

    getCenter(){
        let n = this.vertices.length;
        let [cx,cy] = math.mean(this.vertices,0);
        return [cx,cy];
    }

    getEdges(){
        let n = this.vertices.length;
        let edges = [];
        for (let i = 0;i<n;i++){
            let edge = [this.vertices[i],this.vertices[(i+1)%n]];
            edges.push(edge);
        }
        return edges;
    }

    getEdgeMidpoint(edge_index){
        let edge = this.edges[edge_index];
        let [x1,y1] = edge[0];
        let [x2,y2] = edge[1];
        return [(x1+x2)/2,(y1+y2)/2];
    }

    getNormal(edge_index){
        let edge = this.edges[edge_index];
        let [x1,y1] = edge[0];
        let [x2,y2] = edge[1];
        let [dx,dy] = [x2-x1,y2-y1];
        let normal = [-dy,dx];
        let normal_norm = math.norm(normal);
        normal = [normal[0]/normal_norm,normal[1]/normal_norm];
        //we normalize the normal
        return normal
    }

    getNormals(){
        let n = this.edges.length;
        let normals = [];
        for (let i = 0;i<n;i++){
            let edge_index = i;
            let normal = this.getNormal(edge_index);
            normals.push(normal);
        }

        return normals;
    }



    getClosestPoint2Edge(point,edge){
        //point is an array [x,y]
        //edge is an array of two points [[x1,y1],[x2,y2]]
        //checks if the point projected     
        //on the edge is inside the edge
        //if not returns the closest point on the edge
        let [xp1,xp2] = [point[0]-edge[0][0],point[1]-edge[0][1]];
        let tvec = math.subtract(edge[1],edge[0]);
        let edge_l = math.norm(tvec);
        tvec = math.divide(tvec,edge_l);
        let projection = math.dot(tvec,[xp1,xp2]);
        if (projection<0){
            return [edge[0][0],edge[0][1]];
        }
        if (projection>edge_l){
            return [edge[1][0],edge[1][1]];
        }
        else{
            return math.add(math.multiply(projection,tvec),edge[0]);
        }
    }

    getClosestEdgeNormalProjection(point){
        //Returns \vec{n}^{T}\vec{r_{ep}} where \r_{ep} is the vector from an edge point to the particle (we'll use the midpoint)
        // This is basically the signed distance from the edge to the particle
        // positive if the partice is outside the polygon
        let [closest_edge,distance] = this.getClosestEdge(point);
        let midpoint = this.getEdgeMidpoint(closest_edge);
        let r_ep = math.subtract(point,midpoint);
        let normal = this.getNormal(closest_edge);
        let dot_n_r = math.dot(normal,r_ep);
        return dot_n_r;
    }

    getClosestEdgeNormalProjectionVector(point){
        //returns \vec{n}\vec{n}^{T}\vec{r_{ep}} where \r_{ep} is the vector from an edge point to the particle (we'll use the midpoint)
        let dot_n_r = this.getClosestEdgeNormalProjection(point);
        let normalProjectionVector = math.multiply(dot_n_r,normal);
        return normalProjectionVector.toArray();
    }

    getDistance2Edge(point,edge){
        let closest_point = this.getClosestPoint2Edge(point,edge);
        return math.norm(math.subtract(point,closest_point));
    }


    getClosestEdge(point){
        let n = this.edges.length;
        let min_distance = Infinity;
        let closest_edge = null;

        for (let i = 0;i<n;i++){
            let edge = this.edges[i];
            let closest_point = this.getClosestPoint2Edge(point,edge);
            let distance = math.norm(math.subtract(point,closest_point));
            if (distance<min_distance){
                min_distance = distance;
                closest_edge = i;
            }
        }
        return [closest_edge, min_distance];
    }

    isInteriorPoint(point){
        //checks if a point is inside the polygon
        //we'll use the ray casting algorithm
        //we'll cast a ray from the point to the right (y = point[1]) and count the number of intersections
        //if the number of intersections is odd the point is inside the polygon
        //if the edge is vertical ,dx = 0 so we use any of the vertices x as intersection point
        let n = this.edges.length;
        let intersection = 0;
        let x = point[0];
        let y = point[1];

        for (let i = 0;i<n;i++){
            let edge = this.edges[i];
            let [x1,y1] = edge[0];
            let [x2,y2] = edge[1];

            if (y<Math.min(y1,y2) || y>=Math.max(y1,y2)){
                continue;
            }
            let xIntersection = x1 + (y-y1)*(x2-x1)/(y2-y1);
            //because the ray is casted to the right intersection point should be greater than x

            if (x<=xIntersection){
                intersection++;
            }


        }

        return intersection%2 == 1;
    }
        

    }



    class CollisionHandler{
        /*
        *   This class will handle the collisions between particles and polygons
        *   We'll store the collisions in a map with key value pairs 
        *  iP_jE_Pid:[iP,jE,polygon]  where iP,jE,Pid are the particle index, edge index and polygon id
        * and polygon is the polygon object with id Pid
        *  Methods:
        * detectNewCollisions(collisions, xs, polygons, treshold) -> returns a map with the new collisions
        * detectResolvedCollisions(collisions, xs, treshold) -> returns a map with the resolved collisions
        * updateCollisionsMap(collisions, xs, polygons, treshold) -> updates the collisions map, main function
        * Intended use example in the main loop:
        * let collisions = new Map();
        * STATE->xs, polygons, treshold collisions (besides other properties)
        * let collisionHandler = new CollisionHandler();
        * 
        * let [collisions,new_collisions,resolved_collisions] = collisionHandler.updateCollisionsMap(collisions, xs, polygons, treshold);
        * 
        * 
        * 
        * 
        * 
        * 
        */
        constructor(){
            // We'll store collisions as key value pairs, with id particle_index_edge_index_polygon_id
            // this way it will be easy to check if a collision is already in the map and to remove if
            // it's been resolved

        }
        _detectNewCollisions(collisions, xs, polygons, threshold){
            //brute force method to detect new collisions
            let n = xs.length;
            let m = polygons.length;    


            //the colliding particles
            let new_collisions = new Map();
           
            for (let i = 0;i<n;i++){
                for (let j = 0;j<m;j++){
                    
                    let point = xs[i];
                    let polygon = polygons[j];
                    let [closest_edge,distance] = polygon.getClosestEdge(point);

                    //if the distance is less than the threshold we'll add the collision
                    //additionally to prevent tunneling we check if the particle is inside the polygon
                    // so the conditions to add collision is isInside or distance<threshold
                    let condition = distance<threshold || polygon.isInteriorPoint(point);
                    if (condition){
                        let collision_id = `${i}_${closest_edge}_${polygon.id}`;
                        //if the collision is not already in the map we'll add it
                        if (!collisions.has(collision_id)){        
                            new_collisions.set(collision_id,[i,closest_edge,polygon]);
                            console.log(`NEW COLLISION :Particle ${i} with position ${point} closest edge ${closest_edge} distance ${distance} threshold ${threshold}`);

                        };
                    };
                        }
            
                    }
            return new_collisions;
                };
            


        _detectResolvedCollisions(collisions, xs , treshold){
            //we'll give the resolved collisions, the ones that are not in contact anymore

            let resolved_collisions = new Map();
            let collisions_keys = Array.from(collisions.keys());
            for (let i = 0;i<collisions_keys.length;i++){
                let collision = collisions.get(collisions_keys[i])
                let [particle_index,edge_index,polygon] = collision;
                let collision_id = `${particle_index}_${edge_index}_${polygon.id}`;
                
                // we check if the particle is inside the polygon
                // if it is we dont remove the collision
                // if it is not in the interior we find the closest distance to the edge
                // if the distance is greater than the threshold we remove the collision
                let point = xs[particle_index];
                if (!polygon.isInteriorPoint(point)){
                    let current_distance = polygon.getDistance2Edge(point,polygon.edges[edge_index]);
                    
                    if (current_distance>treshold){
                        resolved_collisions.set(collision_id,collision);
                        console.log(`RESOLVED COLLISION:Particle ${particle_index} with position ${xs[particle_index]} edge ${edge_index} distance ${current_distance} treehold ${treshold}`);
                    }
                }
            }

           return resolved_collisions;   
        }
    

        _updateCollisionsMap(collisions, xs, polygons, threshold){
            //collisions is a map with the current collisions with key [particle_index_edge_index_polygon_id]
            
            let new_collisions = this._detectNewCollisions(collisions, xs, polygons, threshold);
            let resolved_collisions = this._detectResolvedCollisions(collisions, xs, threshold);
            let new_collisions_keys = Array.from(new_collisions.keys());

            for (let i = 0;i<new_collisions_keys.length;i++){
                // lets throw an exception if the collision is already in the map
                if (collisions.has(new_collisions_keys[i])){
                    throw `Collision ${new_collisions_keys[i]} already in the map`;
                }

                //lets avoid corner issues by allowing only one collision per particle
                //if the particle is already in a collision we'll remove the previous collision
                let [particle_index,edge_index,polygon] = new_collisions.get(new_collisions_keys[i]);
                let keys = Array.from(collisions.keys());
                for (let j = 0;j<keys.length;j++){
                    let [p_index,_1,_2] = collisions.get(keys[j]);
                    if (p_index == particle_index){
                        //lets add the previous collision to the resolved collisions
                        resolved_collisions.set(keys[j],collisions.get(keys[j]));
                        collisions.delete(keys[j]);
                    }
                }
                collisions.set(new_collisions_keys[i],new_collisions.get(new_collisions_keys[i]));
            }


            let resolved_collisions_keys = Array.from(resolved_collisions.keys());
            for (let i = 0;i<resolved_collisions_keys.length;i++){
                
                collisions.delete(resolved_collisions_keys[i]);
            }
            
            return [collisions,new_collisions,resolved_collisions];
         
        }

        updateCollisions(STATE, threshold){
            //this method will be the entry point for the collision detection
            //we just have to adapt the state collisionArray and polygonVertices to the required format
            //collisionArray in state is an array of [particle_index, edge_index, polygon_idx]
            //polygonsVertices is an array of vertices arrays



            let [_1,_2,polygonObjs] = stepUnpackStateWithContact(STATE);
            let collisionsArray = STATE.collisions;
            let xs = STATE.xs;
            let collisionsMap = new Map();
            collisionsArray.forEach((collision) => {
                let key = collision[0]+"_"+collision[1]+"_"+collision[2];
                let value = [collision[0],collision[1],polygonObjs[collision[2]]];
                collisionsMap.set(key,value);
            });          

            let [colMap, newcolMap,rescolMap] = this._updateCollisionsMap(collisionsMap, xs, polygonObjs, threshold);

            let collArray = Array.from(colMap.values());
            let newCollArray = Array.from(newcolMap.values());
            let resolvedCollArray = Array.from(rescolMap.values());
            //last element of the arrays is a polygon object, we'll replace it with the polygon id
            collArray = collArray.map((coll) => [coll[0],coll[1],coll[2].id]);
            newCollArray = newCollArray.map((coll) => [coll[0],coll[1],coll[2].id]);
            resolvedCollArray = resolvedCollArray.map((coll) => [coll[0],coll[1],coll[2].id]);



            return [collArray, newCollArray, resolvedCollArray];
        }


    }



//Lets implement the particle-polygon contact constraint,
//This will be modularized in constraints.js but we'll implement it here for now
// the class is intented to standarize the treatment of contact constraints as the other constraints

class ConstraintContact{
    // this class assumes that the constraint is active
    // is the outer loop responsibility to check if the constraint is active and remove the resolved constraints
    //so the edge won't change while the constraint is active by design
    constructor(particle_index, edge_index, polygon){
        this.particle_index = particle_index;
        this.edge_index = edge_index;
        this.polygon = polygon;
    }

    getType(){
        return "contact";
    }

    getParticleIndices(){
        return [this.particle_index];
    }

    getConstraintValue(xarray){
        //C(x) = 0 if the constraint is satisfied
        //It is basically the signed distance from the particle to the edge \vec{n}^{T}\vec{r_{ep}}
        let point = xarray[this.particle_index];
        let dot_n_r = this.polygon.getClosestEdgeNormalProjection(point);
        return dot_n_r;

    }

    getJacobianParticles(xarray){
        //Jacobian is the normal vector of the edge
        let normal = this.polygon.getNormal(this.edge_index);
        return [normal];
    }

    getDotJacobianParticles(varray){
        //dot_Jacobian is zero (static polygon, piecewise linear)

        return [[0,0]];
    }
}

//lets implement computeJacobians function, this will be an extension of the previous computeJacobiansNonContact
// but we'll deal with distance, pin and contact constraints
//order of constraints is [distance, pin, contact]
//if future constraints are added, contact should be the last one
//because they are removed/added as collisions are resolved/detected

function computeJacobians(constraints,xarray,varray){
    /*
    *   Compute the jacobian matrix J and dot jacobian matrix dot_J for the constraints
    *   J is a matrix of size [nconstraints, 2*nparticles]
    *  dot_J is a matrix of size [nconstraints, 2*nparticles]
    *  The constraints should be an array of constraint objects
    */

    // type checking, constraints,xarray,varray should be arrays
    if (!Array.isArray(constraints) || !Array.isArray(xarray) || !Array.isArray(varray)){
        throw new Error("constraints, xarray and varray should be arrays")
    }
    const nconstraints = constraints.length;
    const nparticles = xarray.length;
    const J = math.zeros(nconstraints, 2*nparticles,"sparse");
    const dot_J = math.zeros(nconstraints, 2*nparticles,"sparse");

    for (let i = 0; i < nconstraints; i++){
        let constraint = constraints[i];
        let constraint_particle_indices = constraint.getParticleIndices();
        let jacobian_particles = constraint.getJacobianParticles(xarray);
        let dot_jacobian_particles = constraint.getDotJacobianParticles(varray);

        if (constraint_particle_indices.length !== jacobian_particles.length){
        throw new Error("The number of particles in the constraint is not equal to the number of jacobian particles")
        }
        if (constraint_particle_indices.length !== dot_jacobian_particles.length){
        throw new Error("The number of particles in the constraint is not equal to the number of dot jacobian particles")
        }


        for (let j = 0; j < constraint_particle_indices.length; j++){
            let particle_index = constraint_particle_indices[j];
            let jacobian_particle = jacobian_particles[j];
            if (jacobian_particle.length !== 2){
                throw new Error("The jacobian particle should have 2 components")
            }
            let jacobian_particle_x = jacobian_particle[0];
            let jacobian_particle_y = jacobian_particle[1];
            J.set([i, 2*particle_index], jacobian_particle_x);
            J.set([i, 2*particle_index + 1], jacobian_particle_y);

            let dot_jacobian_particle = dot_jacobian_particles[j];
            if (dot_jacobian_particle.length !== 2){
                throw new Error("The dot jacobian particle should have 2 components")
            }
            let dot_jacobian_particle_x = dot_jacobian_particle[0];
            let dot_jacobian_particle_y = dot_jacobian_particle[1];
            dot_J.set([i, 2*particle_index], dot_jacobian_particle_x);
            dot_J.set([i, 2*particle_index + 1], dot_jacobian_particle_y);

        
        }
    }

    return [J, dot_J];
    }





//lets implement the contact constraint solver, this will be an extension of the previous solveNonContactConstraints
// In principle there is nothing to change, it is a copy paste, contact constraints are treated as any other constraint
// it is the outer step loop that will have some changes, like clipping the contact lagrange multipliers 
//(they can't be negative, contact forces are repulsive) and removing the contact constraints when the collision is resolved
function calculateConstraintForces(constraints, lambda_multipliers, J){
    //this function will split the constraint forces in the different types of constraints
    //mainly for debugging purposes
    //every constraint has a getType method
    //we separate the constraint forces by type
    //constraints is a [nconstraints] of constraint objects
    //lambda_multipliers is [nconstraints] array of lagrange multipliers
    // J is a [nconstraints,2*nparticles] matrix
    // so we iterate over the constraints, get the lambda_constr_type_I and the J_constr_type_I for all constraints
    // and we return the J_constr_type_I^{T}*lambda_constr_type_I 
    // additionally we sum all the constraint forces and return it as well
    // so we must return an object with keys [distance,pin,contact,all]

    let nconstraints = constraints.length;
    let nparticles = J.size()[1]/2;
    let constraintForces = math.zeros(nparticles,2);
    let constraintForcesDistance = math.zeros(nparticles,2);
    let constraintForcesPin = math.zeros(nparticles,2);
    let constraintForcesContact = math.zeros(nparticles,2);

    for (let i = 0; i < nconstraints; i++){
        let constraint = constraints[i];
        let lambda = lambda_multipliers[i];
        let J_i = J.subset(math.index(i,math.range(0,2*nparticles)));
        let constraintForce = math.multiply(math.transpose(J_i), lambda).reshape([nparticles,2]);

        constraintForces = math.add(constraintForces, constraintForce);
        if (constraint.getType() === "distance"){
            constraintForcesDistance = math.add(constraintForcesDistance, constraintForce);
        }
        if (constraint.getType() === "pin"){
            constraintForcesPin = math.add(constraintForcesPin, constraintForce);
        }
        if (constraint.getType() === "contact"){
            constraintForcesContact = math.add(constraintForcesContact, constraintForce);
        }
    }

    return {
        distance_forces: constraintForcesDistance.toArray(),
        pin_forces: constraintForcesPin.toArray(),
        contact_forces: constraintForcesContact.toArray(),
        constraint_forces: constraintForces.toArray(),
    }
}
function solveConstraints(J, dot_J, xarray, varray,
    external_forces, masses, constraints,
     alpha = 0.00, beta = 0.00){


    //lets do type checking, J,dot_J external forces mathjs matrices, xarray, varray, masses arrays
    if (!math.isMatrix(J) || !math.isMatrix(dot_J) || !math.isMatrix(external_forces)){
    throw new Error("J, dot_J and external_forces should be mathjs matrices")
    }
    if (!Array.isArray(xarray) || !Array.isArray(varray) || !Array.isArray(masses)){
    throw new Error("xarray, varray and masses should be arrays")
    }
    let nparticles = xarray.length;
    let nconstraints = J.size()[0];
    let M_inv = math.matrix(masses.map(mass=>[1/mass,1/mass])).reshape([-1]);//vx and vy have the same mass(same particle)
    M_inv = math.diag(M_inv,"sparse");

    let A = math.multiply(J, math.multiply(M_inv, math.transpose(J)));
    //lets check shape of A, it should be [nconstraints,nconstraints]
    if (A.size()[0] !== nconstraints || A.size()[1] !== nconstraints){
    throw new Error("The A matrix should have the same number of constraints")
    }
    //lets compute b = -dot_J*varray - J*M_inv*external_forces - alpha*J*varray - beta*constraints
    //Lets do it in parts and then sum
    //first varray [nparticles,2] to matrix and flatten-> [2*nparticles]
    //the same for external forces [nparticles,2] to matrix and flatten-> [2*nparticles]
    //Get constraints_vals [nconstraints] from constraints

    let constraints_vals = math.matrix(
        constraints.map(constraint => constraint.getConstraintValue(xarray))
        );
    let vmatrix = math.flatten(varray);
    external_forces =math.flatten(external_forces);
    let part1 = math.multiply(dot_J, vmatrix);
    let part2 = math.multiply(J, math.multiply(M_inv, external_forces));
    let part3 = math.multiply(alpha, math.multiply(J, vmatrix));
    let part4 = math.multiply(beta, constraints_vals).reshape([nconstraints,1]);


    let b = math.add(math.add(part1, part2), math.add(part3, part4));
    //all parts appear with a minus sign in the equation
    b = math.multiply(b,-1);
    //lets add for stability a small value to the diagonal of A
    //let I = math.identity(nconstraints);
    //A = math.add(A, math.multiply(1e-5,I));

    let lagrange_multipliers = math.lusolve(A,b).reshape([-1]);
    // We have to clip the contact lagrange multipliers to be >= 0
    // lets get the indices of the contact constraints using the getType method
    let contact_indices = constraints.map((constraint,i) => [constraint.getType(),i]).filter(c => c[0] === "contact").map(c => c[1]);
    for (let i = 0; i < contact_indices.length; i++){
        let index = contact_indices[i];
        if (lagrange_multipliers.get([index]) < 0){
            lagrange_multipliers.set([index],0);
        }
    }

    //let constraintForces = math.multiply(math.transpose(J), lagrange_multipliers).reshape([nparticles,2]).toArray();

    let forcesDict = calculateConstraintForces(constraints, lagrange_multipliers.toArray(), J)

    let out_dict = {
        constraintForces: forcesDict.constraint_forces,
        lagrange_multipliers: lagrange_multipliers.toArray(),
        constraints_vals: constraints_vals.toArray(),
        contact_forces: forcesDict.contact_forces,
        pin_forces: forcesDict.pin_forces,
        distance_forces: forcesDict.distance_forces,

    }
    return out_dict;
     }

//lets write a util function to transform the collision map to a coontact_constraints object array
//In the State object, collisions will be an array of [iP,jE,Pid] where iP is the particle index, jE is the edge index and Pid is the polygon id
// We'll use polygon index in the polygons array as polygon id , so collisions will be an array of [iP,jE,Pidx]
// the polygons array will be an array of vertices arrays, so we'll create the Polygon objects here

function collisions2ContactConstraints(collisions, polygonsVertices){
    let contactConstraints = [];
    for (let [iP,jE,Pidx] of collisions){
        let polygon = new Polygon(polygonsVertices[Pidx],Pidx);
        let contactConstraint = new ConstraintContact(iP,jE,polygon);
        contactConstraints.push(contactConstraint);
    }
    return contactConstraints;
}

function calculateFrictionForces(lagrangeMultipliersContact, varray, collisions, polygonsObjs,muFriction){
    // friction force is = -mu*normal_force*Vt/|Vt|
    //V_{particleTangent} = V_{particle} - \vec{n}^{T}V_{particle}\vec{n} with \vec{n} the normal vector of the edge, we get it from the polygon object
    //collisions is a [[iP,jE,Pidx],...] array where iP is the particle index, jE is the edge index and Pidx is the polygon index
    // There must be the same lagrangeMultipliersContact as the number of collisions
    
    let frictionForces = math.zeros(varray.length,2);
    //throw an error if the number of lagrange multipliers is different from the number of collisions
    if (lagrangeMultipliersContact.length !== collisions.length){
        throw new Error("The number of lagrange multipliers should be equal to the number of collisions")
    }

    for (let i = 0; i < collisions.length; i++){
        let [iP,jE,Pidx] = collisions[i];
        let polygon = polygonsObjs[Pidx];
        let normal = polygon.getNormal(jE);
        let Vparticle = varray[iP];
        let VparticleTangent = math.subtract(Vparticle, math.multiply(math.dot(normal,Vparticle),normal));
        let Vt = math.norm(VparticleTangent);
        let normalForce = lagrangeMultipliersContact[i];
        let eps = 1e-6;
        let frictionForce = math.multiply(-muFriction*normalForce/(Vt+eps),VparticleTangent);
        frictionForces.set([iP,0], frictionForce[0]);
        frictionForces.set([iP,1], frictionForce[1]);

    
    }

    return frictionForces.toArray();
}


// Lets write a function to unpack the state object. This function in the future will replace the stepUnpackState function
// but for now lets write an additional one to avoid breaking the code

function stepUnpackStateWithContact(STATE){
    // We turn the constraint info into constraint objects
    // additionally we turn force info into force objects
    // I'll enumerate the field and format  of the different entries in STATE, relevant for this function
    //constraints_distance: [[idx1, idx2, distance],...] ConstraintDistance class has signature  constructor(particle1_idx, particle2_idx, distance)
    // constraints_pin: [[idx,pin_point, pin_distance],...] ConstraintPin class has signature constructor(particle_idx, pin_point, pin_distance)
    // collisions: [[iP,jE,Pidx],...] where iP is the particle index, jE is the edge index and Pidx is the polygon index
    // polygons: [[vertex1,vertex2,...],...]
    // constant_forces: [[idx1, force1],...] 
    // gravity: [g, direction] Gravity class has signature constructor(particle_indices, g, direction)
    // by default particles_indices arreys are [0,1,2,...,Nparticles-1]
    // damping: [[idx1, damping_coefficient],...] Damping class has signature constructor(particle_indices, damping_coefficient)
    // springs: [[idx1,idx2, spring_constant, rest_length],...] Spring class has signature constructor(idx1,idx2, spring_constant, rest_length)
    // springs2points: [[idx1,fixed_point, spring_constant, rest_length],...] Spring2Point class has signature constructor(idx1,fixed_point, spring_constant, rest_length) 

    //lets throw an error if the following keys are not present
    let mandatory = ["constraints_distance", "constraints_pin","collisions","polygons", "xs", "vs", "masses"]
    for (let key of mandatory){
        if (!STATE.hasOwnProperty(key)){
            throw new Error(`The STATE object must have the key ${key}`)
        }
    }

    let s = STATE;

    let constraintsDistanceObjs = s.constraints_distance.map(c => new ConstraintDistance(c[0], c[1], c[2]));
    let constraintsPinObjs = s.constraints_pin.map(c => new ConstraintPin(c[0], c[1], c[2]));
    let contactConstraints = collisions2ContactConstraints(s.collisions, s.polygons);
    let polygonsObjs = s.polygons.map((p,i) => new Polygon(p,i));
    let constraintsObjs = constraintsDistanceObjs.concat(constraintsPinObjs).concat(contactConstraints);
    let nparticles = s.xs.length;
    let nparticles_array = Array.from(Array(nparticles).keys());
    let forceObjs = [];
    //if constant_forces/gravity/damping/springs not present, we dont add nothing to forceObjs
    if (s.hasOwnProperty("constant_forces")){
        forceObjs = s.constant_forces.map(f => new ConstantForce([f[0]], f[1]));
    }
    if (s.hasOwnProperty("gravity")){
        forceObjs.push(new Gravity(nparticles_array, s.gravity[0], s.gravity[1]));
    }
    if (s.hasOwnProperty("damping")){
        forceObjs = forceObjs.concat(s.damping.map(f => new Damping([f[0]], f[1])));
    }
    if (s.hasOwnProperty("springs")){
        forceObjs = forceObjs.concat(s.springs.map(f => new Spring(f[0], f[1], f[2], f[3])));
    }
    if (s.hasOwnProperty("springs2points")){
        forceObjs = forceObjs.concat(s.springs2points.map(f => new Spring2Point(f[0], f[1], f[2], f[3])));
    }



    return [constraintsObjs, forceObjs, polygonsObjs];

}
    


//lets implement the step function with contact constraints

function frictionForcesTweak(frictionForces,J){
    // we'll tweak the friction forces to be more consistent with the constraints
    // basically we want to sustract to teh friction forces its projection on the constraints normal vectors
    // this will prevent the friction forces from pushing the particle into the polygon
    // we want Jf = 0, where f is the friction forces, to avoid solving the system again we'll just tweak the friction forces
    let frictionForcesFlat = math.matrix(math.flatten(frictionForces));
    // we iter over each row of J^{T} and substract the projection of the friction force on the normal vector
    let nconstraints = J.size()[0];
    let nparticles = J.size()[1]/2;
    let frictionForcesTweaked = math.clone(frictionForcesFlat);
    let Jarray = J.toArray();
    for (let i = 0; i < nconstraints; i++){
        let Jt_i = Jarray[i];

        let normJt_i = math.norm(Jt_i);
        let projection = math.dot(Jt_i, frictionForcesFlat)/normJt_i;
        frictionForcesTweaked = math.subtract(frictionForcesTweaked, math.multiply(projection,Jt_i));

    }

    return frictionForcesTweaked.reshape([nparticles,2]);
}



function stepSemiEulerContact(STATE, PARAMETERS){
    // things we have to do here:
    // 1. Unpack the state object
    // 2. Compute the Jacobians 
    // 3. Calculate the external forces
    // 4. Solve the constraints
    // 5. We calculate the friction force once we know the contact normal force
    // 6. Integrate the system
    // 7. Update collisions

    let mandatoryParameters = ["collisionThreshold", "muFriction", "alpha", "beta"];
    for (let key of mandatoryParameters){
        if (!PARAMETERS.hasOwnProperty(key)){
            throw new Error(`The PARAMETERS object must have the key ${key}`)
        }
    }

    let [constraintObjs, forceObjs, polygonsObjs] = stepUnpackStateWithContact(STATE);
    let s = STATE;
    let [J,dotJ] = computeJacobians(constraintObjs, s.xs, s.vs);
    let externalForces = computeExternalForces(forceObjs, s.xs, s.vs, s.masses);
    let outSolverConstraints = solveConstraints(J, dotJ, s.xs, s.vs, externalForces, s.masses, constraintObjs, PARAMETERS.alpha, PARAMETERS.beta);
    let indicesConstraintsContact = constraintObjs.map((c,i) => [c.getType(),i]).filter(c => c[0] === "contact").map(c => c[1]);
    let lagrangeMultipliersContact = indicesConstraintsContact.map(i => outSolverConstraints.lagrange_multipliers[i]);

    let frictionForces = calculateFrictionForces(lagrangeMultipliersContact, s.vs, s.collisions, polygonsObjs, PARAMETERS.muFriction);
   // frictionForces = frictionForcesTweak(frictionForces,J).toArray();
    let constraintForces = outSolverConstraints.constraintForces;
    let contactForces = outSolverConstraints.contact_forces;
    let totalForces = math.add(math.add(externalForces,constraintForces), frictionForces).toArray();
    let acc  = totalForces.map((f,i) => [f[0]/s.masses[i], f[1]/s.masses[i]]);


    let dt = PARAMETERS.dt;
    let newVs = math.add(s.vs, math.multiply(dt,acc));
    let newxs = math.add(s.xs, math.multiply(dt,newVs));

    let colHandler = new CollisionHandler();

    //newState with the values updated in the step
    let newState = {
        "xs": newxs,
        "vs": newVs,
        "time": s.time + dt,
        "friction_forces": frictionForces,
        "external_forces": externalForces.toArray(),
        "constraint_forces": constraintForces,
        "contact_forces": contactForces,
        "total_forces": totalForces,
        "lagrange_multipliers": outSolverConstraints.lagrange_multipliers,
        "constraints_vals": outSolverConstraints.constraints_vals,
        "J": J.toArray(),
        "dot_J": dotJ.toArray(),

    }
    let [collisionsArray, newCollisionsArray, resolvedCollisionsArray] = colHandler.updateCollisions(STATE, PARAMETERS.collisionThreshold);
    
    newState.collisions = collisionsArray;
    newState.new_collisions = newCollisionsArray;
    newState.resolved_collisions = resolvedCollisionsArray;


    //lets add to newState the state properties not included in the previous object
    //we loop over STATE keys and if they are not present in newState we add them
    for (let key in STATE){
        if (!newState.hasOwnProperty(key)){
            newState[key] = STATE[key];
        }
    }

    return newState;
}

























/*
TESTCASES
*/
// Here we'll specify some simple systems with analytical solutions to test the solver
// we have to specify the initial STATE object, and the analytical solution function
// the analytical solution function should be a function of x0, v0, t 
// we'll start with a Pendulum and a Dumbell




function calculateConstraintForcesFromTrajectory(vsArray, masses, externalForces,tArray){
    //utility to numerically calculate constraint forces from the trajectory velocities calculated by the solver or the analytical solution
    // constraint_forces are nothing but Mx''-externalForces = constraint_forces, with M the mass matrix, a diagonal matrix in our case
    //varray is a [tsteps, nparticles, 2] array with the velocities masses [nparticles] externalForces [tsteps,nparticles,2]
    //we assume fixed size dt
    let shapesEqual = isShapeEqual(vsArray, externalForces);
    if (!shapesEqual){
        throw new Error("The shapes of vArray and externalForces should be equal")
    }
    
    let h = tArray[1]-tArray[0];

    let accArray = [];
    for (let i = 0; i < vsArray.length-1; i++){
        let acc = vsArray[i+1].map((v,j) => [(v[0]-vsArray[i][j][0])/h, (v[1]-vsArray[i][j][1])/h]);
        accArray.push(acc);
    }

    accArray.push(accArray[accArray.length-1]);


    let constraintForces = [];
    for (let i = 0; i < tArray.length; i++){
        let acc = accArray[i];
        let externalForce = externalForces[i];
        let constraintForce = acc.map((a,j) => [masses[j]*a[0]-externalForce[j][0], masses[j]*a[1]-externalForce[j][1]]);
        constraintForces.push(constraintForce);
    }

    return constraintForces;

}






function pendulumSystemConfig(thetaO, omegaO, g, L, mass){

    let x0 = [[L*Math.sin(thetaO), -L*Math.cos(thetaO)]]; // [nparticles,2] [1,2] ( 1 particle)
    let v0 = [[omegaO*L*Math.cos(thetaO), omegaO*L*Math.sin(thetaO)]];
    let masses = [mass];
    let constraints_distance = [];
    let constraints_pin = [[0, [0,0], L]];
    let gravity = [g, [0,-1]];
    let time = 0;
    let initialState = {
        xs: x0,
        vs: v0,
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: constraints_pin,
        collisions: [],
        gravity: gravity,
        time: time,
        polygons: [],
        }

    
    
    return initialState;

}

function pendulumAnalyticalSolution(thetaO, g, L, tarray){
    //small angle approximation
    // lets return a [tarray.length,2] array with the positions of the particle at each time
    // the same for vs
    let omega = Math.sqrt(g/L);
    let theta = thetaO*Math.cos(omega*tarray);
    let omega_t = -thetaO*omega*Math.sin(omega*tarray);
    let x = L*Math.sin(theta);
    let y = -L*Math.cos(theta);
    let vx = omega*L*Math.cos(theta);
    let vy = omega*L*Math.sin(theta);
    let xs = [];
    let vs = [];
    for (let t of tarray){
        let theta = thetaO*Math.cos(omega*t);
        let omega_t = -thetaO*omega*Math.sin(omega*t);
        let x = L*Math.sin(theta);
        let y = -L*Math.cos(theta);
        let vx = omega*L*Math.cos(theta);
        let vy = omega*L*Math.sin(theta);
        xs.push([x,y]);
        vs.push([vx,vy]);
    }
    return [xs,vs];
}

function doublePendulumSystemConfig(thetaO1, omegaO1, thetaO2, omegaO2, g, L1, L2, mass1, mass2){
    let x00 = [L1*Math.sin(thetaO1), -L1*Math.cos(thetaO1)];
    let x01 = [x00[0] + L2*Math.sin(thetaO2), x00[1] - L2*Math.cos(thetaO2)];
    let x0 = [x00, x01];
    let v0 = [[omegaO1*L1*Math.cos(thetaO1), omegaO1*L1*Math.sin(thetaO1)], [omegaO2*L2*Math.cos(thetaO2), omegaO2*L2*Math.sin(thetaO2)]];
    let masses = [mass1, mass2];
    let constraints_distance = [ [0,1, L2]];
    let constraints_pin = [[0, [0,0], L1]];
    let gravity = [g, [0,-1]];
    let time = 0;
    let initialState = {
        xs: x0,
        vs: v0,
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: constraints_pin,
        gravity: gravity,
        time: time,
        }

    return initialState;
}

function chainSystemConfig(nparticles, g, L, mass, angle){
    // we'll create a chain of particles with a given angle to avoid an initial static configuration
    let x0 = [];
    let v0 = [];
    let masses = [];
    let constraints_distance = [];
    let constraints_pin = [];
    let gravity = [g, [0,-1]];
    let time = 0;
    constraints_pin.push([0, [0,0], L]);
    for (let i = 0; i < nparticles; i++){

        // lets add a little angle perturbation to avoid singular configurations
        let angle_perturbation = 0.01*Math.random()+angle;
        let x = [(i+1)*L*Math.cos(angle_perturbation), -(i+1)*L*Math.sin(angle_perturbation)];
        let v = [0,0];
        x0.push(x);
        v0.push(v);
        masses.push(mass);
        if (i > 0){
            constraints_distance.push([i-1,i,L]);
        }
    }

    let initialState = {
        xs: x0,
        vs: v0,
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: constraints_pin,
        gravity: gravity,
        time: time,
        }

    return initialState;
}


// With contact systems

function particleOnStairs(widthStep =0.3, heightStep = 0.05,vx = 0.1,nsteps =3, gravity = 0.1, mass = 1){
    //a simple scenario to test contact constraints and resolution
    //n steps, the particle initiall rests on the higher one, but moving
    //with a given tangential velocity, it must move with constant velocity
    // until it reaches the end of the step, for per4flectly inelastic collisions
    //it must lost all the normal velocity and keep the tangential velocity

    //for the stairs we just create rectangles polygons widthStep , heighStep
    //and we shift them widthstep, heightstep
    // the particle will be at the center of the first step

    let polygons = [];
    let heightPol = 0.1;
    for (let i = 0; i < nsteps; i++){
        let tLC = [i*widthStep,- i*heightStep];//top left corner
        let polygon = [tLC, [tLC[0]+widthStep, tLC[1]],
                     [tLC[0]+widthStep, tLC[1]-heightStep], [tLC[0], tLC[1]-heightStep] ];

        polygons.push(polygon);
    }

    let x0 = [[widthStep/2, 0]];
    
    let initialState = {
        xs: x0,
        vs: [[vx,0]],
        masses: [mass],
        constraints_distance: [],
        constraints_pin: [],
        collisions: [],
        gravity: [gravity, [0,-1]],
        time: 0,
        polygons: polygons,
        }

    return initialState;


}

function dumbellSlidingOnWallSystemConfig(L,mass1,mass2,angle, g){
    //simple test for contact constraints, without friction
    // a dumbell system, two masses with a distance constraint, mass1 is on the ground
    //mass2 is on the wall, angle is the initial angle of the dumbell with the horizontal
    //the system will slide down the wall with an effective 1 degree of freedom (the angle)
    //we must create two polygons, one for the ground and one for the wall
    // and we place mass1 at L(cos angle,0) and mass2 at L(0, sin angle)

    // we build polygons clock-wisely so the edges and normals are correctly computed

    let polWall = [[0,-3],[-1,-3],[-1,3],[0,3]];
    let polGround = [[-3,0],[3,0],[3,-1],[-3,-1]];
    angle = angle*Math.PI/180;
    let x1 = [L*Math.cos(angle),0];
    let x2 = [0,L*Math.sin(angle)];
    let x0 = [x1,x2];
    let v0 = [[0,0],[0,0]];
    let masses = [mass1,mass2];

    let initialState = {
        xs: x0,
        vs: v0,
        masses: masses,
        constraints_distance: [[0,1,L]],
        constraints_pin: [],
        collisions: [],
        gravity: [g,[0,-1]],
        time: 0,
        polygons: [polGround,polWall],
        }



    return initialState;
}

function dumbellSlidingOnWallAnalyticalSolution(L, mass, angle, g, tArray){
    // we'll integrate the angle evolution equation for the simplified assumption of equal masses
    //Lagrangian = T-V. T = 1/2 m1 v1^2 + 1/2 m2 v2^2   and V = mgx2 = mgLsin(theta)
    //v1 = d_t(Lcos(angle)) = -Lsin(angle)theta_dot v2 = d_t(Lsin(angle)) = Lcos(angle)theta_dot
    //T = 1/2 m1 L^2 sin^2(theta) theta_dot^2 + 1/2 m2 L^2 cos^2(theta) theta_dot^2 = mL^{2}theta_dot^2/2
    //V = mgLsin(theta)
    //L = mL^{2}theta_dot^2 - mgLsin(theta), eofm d_t(partial_L/partial_theta_dot) - partial_L/partial_theta = 0
    //2mL^{2}theta_dot theta_ddot + mgLcos(theta) = 0
    //theta_ddot = -g/(2L)cos(theta)

    let theta0 = angle*Math.PI/180;

    let dxdtFun = (x,t) => {
        let [theta,theta_dot] = x;
        let dthetadt = theta_dot;
        let dtheta_dotdt = (-g/(L))*Math.cos(theta);
        return [dthetadt,dtheta_dotdt];
    
    }

    let x0 = [theta0,0];
    let out = integrateOdeOnTarray(dxdtFun, x0, tArray,"Euler")//out [tarray.length,2] 
    let thetaArray = out.map(x => x[0]);
    let theta_dotArray = out.map(x => x[1]);

    let xsArray = thetaArray.map((theta,i) => [ [L*Math.cos(theta),0], [0,L*Math.sin(theta)] ]);
    let vsArray = theta_dotArray.map((theta_dot,i) => [ [-L*Math.sin(thetaArray[i])*theta_dot,0], [0,L*Math.cos(thetaArray[i])*theta_dot] ]);
    let externalForces = tArray.map(t => [ [0,-mass*g],[0,-mass*g] ]);
    let constraintForces = calculateConstraintForcesFromTrajectory(vsArray, [mass,mass], externalForces, tArray);

    let stateStory = [];

    for (let i = 0; i < tArray.length; i++){
        let state = {
            "xs": xsArray[i],
            "vs": vsArray[i],
            "time": tArray[i],
            "constraint_forces": constraintForces[i],
            "external_forces": externalForces[i],
            "gravity": [g,[0,-1]],
            "constraints_distance": [[0,1,L]],
            "polygons": [[[-3,0],[3,0],[3,-1],[-3,-1]],[[0,-3],[-1,-3],[-1,3],[0,3]]],
            "masses": [mass,mass],

        }
        stateStory.push(state);
    }

    return stateStory;
    
}

        
function tumblingBoxSystemConfig(side, g, Lwedge, angle){
    //here we aim to simulate the tumbling box along a wedge
    let origin = [0,0];
    let xs = [[origin[0],origin[1]],
              [origin[0],origin[1]+side],
                [origin[0]+side,origin[1]+side],
                [origin[0]+side,origin[1]]];

    let neighs = getNeighborsDelauney(xs);//loaded from math_utils.js

    let visitedPairs = new Set();
    let constraints_distance = [];
    for (let particle_idx in neighs){
        let neighbors = neighs[particle_idx];
        for (let i = 0; i < neighbors.length; i++){
            let particle_idx2 = neighbors[i];
            let pair = [particle_idx, particle_idx2];
            pair.sort();
            let pairStr = pair.join(",");
            if (!visitedPairs.has(pairStr)){
                visitedPairs.add(pairStr);
                //particle_idx and particle_idx2 to Integer
                particle_idx = parseInt(particle_idx);
                particle_idx2 = parseInt(particle_idx2);
                let distance = math.distance(xs[particle_idx], xs[particle_idx2]);
                constraints_distance.push([particle_idx, particle_idx2, distance]);
            }
        }
    }
    
    

    let base = Lwedge;
    //angle in deg
    let height = base*Math.tan(angle*Math.PI/180);
    let wedge = [[0,0],[0,height],[base,0]];
    //now, lets place the box on the wedge, we'll place the pivot box vertex at intersection of the 90º angle of the wedge bisector and the wedge slope
    // for this we start by shifting the box by [-side,0] to the left so the pivot is at [0,0], we then rotate the box around the pivot by angle
    //and finally we shift the box by the vector slope_midpoint+origin. 
    //we'll use translate and rotate functions from math_utils.js
    xs = translate(xs, [-side,0]);
    xs = rotate(xs,angle,[0,0]);
    //midpoint [0,height]+[base/2,height]/2
    let slopeMidpoint = [base/2,height-height/2];

    xs = translate(xs, slopeMidpoint);

    let STATE = {
        neighbors: neighs,
        xs: xs,
        constraints_distance: constraints_distance,
        constraints_pin: [],
        "constraints_values": [],
        "constraints_mae": 0,
        "targets": [],
        "masses": [1,1,1,1],
        "vs": [[0,0],[0,0],[0,0],[0,0]],
        "collisions": [],
        "polygons": [],
        "J": null,
        "dotJ": null,
        "t":0,
        "gravity": [g,[0,-1]],
        "polygons": [wedge], 
        "springs2points": []


    };

    //lets translate box and wedge a little to the left and down
    STATE.xs = translate(STATE.xs, [-1.5,-1]);
    STATE.polygons[0] = translate(STATE.polygons[0], [-1.5,-1]);
    return STATE;
}


function tumblingBoxAnalyticalSolution(side, g, angle, tarray){
    //Here we integrate the angle coordinate evolution
    //we'll focus on the rotation around the pivot point without considering pivot swaps on 
    //following steps. Angle will be positive for clockwise rotation
    //and theta0 = wedge_angle
    //The pivot point will be the 3 vertex of the box (0-indexed) to be
    //consistent with the box configuration we defined in tumblingBoxSystemConfig
    // the equation will be I*theta'' = sum_i torque_i = sum_i r_i x (-m_i*g) = -m*g*sum_i L_i*sin(theta_i)
    // for the 0 and 2 masses L_i = side, the 1 masses is sqrt(2)*side (diagonal to the pivot)
    // angles for the [0,1,2] masses are [theta, theta+pi/4, theta+pi/2]
    // so I*theta'' = -m*g*(L*sin(theta)+L*sqrt(2)*sin(theta+pi/4)+L*sin(theta+pi/2))
    let origin = [0,0];
    let xs = [[origin[0],origin[1]],
            [origin[0],origin[1]+side],
                [origin[0]+side,origin[1]+side],
                [origin[0]+side,origin[1]]];


    let masses = [1,1,1,1];
    let pivot = xs[3];
    let inertiaMoment = calculateInertiaMoment( xs, masses, pivot);
    let theta0 = angle*Math.PI/180;

    //lets write the differential equation as a 2 dimensional first order system
    // so we can use integrateOde(fun, x0, tf, h, method="RK4", tArrayEval = null) from math_utils.js
    //where fun is dxdt(x,t)->[dx1dt,dx2dt] and x0 is [x0_1,x0_2]  
    let dxdt = (x,t) => {
        let theta = x[0];
        let omega = x[1];
        let L = side;
        let [theta1, theta2, theta3] = [theta, theta-Math.PI/4, theta-Math.PI/2];
        let torque = g*L*(Math.sin(theta1)+Math.sqrt(2)*Math.sin(theta2)+Math.sin(theta3));
        let dotTheta = omega;
        let dotOmega = torque/inertiaMoment;
        return [dotTheta, dotOmega];
    }

    let x0 = [theta0, 0];
    let tf = tarray[tarray.length-1];
    let h = tarray[1]-tarray[0];
    let method = "RK4";
    let tArrayEval = tarray;
    let out = integrateOde(dxdt, x0, tf, h, method, tArrayEval);
    let thetas = out.map(x => x[0]);
    let omegas = out.map(x => x[1]);
    // lets transform thetas and omegas to xs and vs

    let xsArray = [];
    let vsArray = [];
    for (let i = 0; i < thetas.length; i++){
        let theta = thetas[i];
        let omega = omegas[i];
        let x = rotate(xs, theta, pivot,"rad");
        //lets get the relative vectors to the pivot (vertex 3)
        let r30 = math.subtract(x[0],pivot);
        let r31 = math.subtract(x[1],pivot);
        let r32 = math.subtract(x[2],pivot);
        //pivot has [0,0] velocity
        let v0 = [r30[1]*omega, -r30[0]*omega];
        let v1 = [r31[1]*omega, -r31[0]*omega];
        let v2 = [r32[1]*omega, -r32[0]*omega];
        let v3 = [0,0];
        let vs = [v0,v1,v2,v3];


        xsArray.push(x);
        vsArray.push(vs);
    }

    //lets calculate the constraint forces
    // basically mx'' = external_forces + constraint_forces
    // so constraint_forces = mx'' - external_forces = mx'' - m*g
    //lets calculate the acceleration numerically acc_t = (v_{t+1}-v_t)/h
    //for the last value we just copy the previous one
    //vsArray is a [tarray.length,4,2] array
    // we have to build a [tarray.length,4,2] array with acceleration values

    let externalForces = [];
    for (let i = 0; i < tarray.length; i++){
        let f1 = [0,-masses[0]*g];
        let f2 = [0,-masses[1]*g];
        let f3 = [0,-masses[2]*g];
        let f4 = [0,-masses[3]*g];

        externalForces.push([f1,f2,f3,f4]);
    }
    let constraintForces = calculateConstraintForcesFromTrajectory(vsArray, masses, externalForces,tArrayEval)
    //lets assert constraintForces are of the same shape as vsArray
    if (constraintForces.length !== vsArray.length){
        throw new Error("constraintForces and vsArray should have the same length")
    }
 




    //now lets create a stateStory array of state objects that mimic the solver output
    let stateStory = [];
    for (let i = 0; i < tarray.length; i++){
        let state = {
            "xs": xsArray[i],
            "vs": vsArray[i],
            "time": tarray[i],
            "constraint_forces": constraintForces[i],
            "external_forces": externalForces[i],
            "masses": masses,
            "constraints_distance": [[0,1,side],[1,2,side],[2,3,side],[3,0,side],[1,3,Math.sqrt(2)*side]],
            "constraints_pin": [],
            "collisions": [],
            "gravity": [g, [0,-1]],
            "polygons": [],
            "J": null,
            "dotJ": null,
            "inertiaMoment": inertiaMoment,//for debugging purposes
        }
        stateStory.push(state);
    }

    return stateStory;


}

function integrateSystemNonContact(STATE,PARAMETERS, nsteps){    
    let s = STATE;

    let STATE_STORY = [s];

    for (let i = 0; i < nsteps; i++){
        //lets copy the state object, just in case
        let state = JSON.parse(JSON.stringify(s));
        let newState = stepSemiEulerNonContact(state,PARAMETERS);
        STATE_STORY.push(newState);
        s = newState;
    }

    return STATE_STORY;



}

function integrateSystem(STATE,PARAMETERS, nsteps){
    let s = STATE;

    let STATE_STORY = [s];

    for (let i = 0; i < nsteps; i++){
        //lets copy the state object, just in case
        let state = JSON.parse(JSON.stringify(s));
        let newState = stepSemiEulerContact(state,PARAMETERS);
        STATE_STORY.push(newState);
        s = newState;
    }

    return STATE_STORY;
}


let referenceConfigurations = {
    "pendulum": {
        "thetaO": Math.PI/4,
        "omegaO": 0,
        "g": 9.8,
        "L": 1,
        "mass": 1,
    },
}
function solveReferenceSystems(nameSystem,PARAMETERS,nsteps){
    //panelPlotsId is the id of the div where the plots will be displayed
    //we'll use plotly to create the plots, plotly will be loaded in the html
    let config = referenceConfigurations[nameSystem];
    let stateStory = [];
    if (nameSystem === "pendulum"){
       
        
        let initialState = pendulumSystemConfig(config.thetaO, config.omegaO, config.g, config.L, config.mass);
        stateStory = integrateSystem(initialState, PARAMETERS, nsteps);


    }

    //lets return the stateStory

    return stateStory;
}


export { stepTestSemiEuler,stepSemiEulerContact, calculateCOM, calculateKineticEnergy,
        stepSemiEulerNonContact, solveReferenceSystems, pendulumSystemConfig, pendulumAnalyticalSolution,
        doublePendulumSystemConfig, integrateSystem, particleOnStairs,
         initState, chainSystemConfig,  tumblingBoxSystemConfig,tumblingBoxAnalyticalSolution,
            dumbellSlidingOnWallSystemConfig, dumbellSlidingOnWallAnalyticalSolution,
        computeJacobians, solveConstraints, stepUnpackStateWithContact,
         collisions2ContactConstraints,computeExternalForces, calculateConstraintForces,
         CollisionHandler, ConstraintContact, Polygon
        };