import { ConstraintDistance, ConstraintPin, ConstraintContact } from "./constraints.js";
import {getAngle, translate, rotate,rotateCOM, getNeighborsDelauney, calculateInertiaMoment,
    integrateOde, integrateOdeOnTarray, isShapeEqual, calculateCOM, calculateKineticEnergy, getNeighborsDict} from "./math_utils.js";
import {ConstantForce,Gravity,Damping,Spring,Spring2Point,computeExternalForces} from "./forces_actuators.js";
//lets import external dependencies: mathjs
import * as math from 'mathjs'; 

//IntegrateSystem(fun, x0, tf, h) where x0 is a R dim array, fun is a R->R dxdt fun 
// integrateOdeOnTarray(fun,xo,tArray,method = "RK4")
//h is the time step, tf is the final time, by default uses RK4
//all constraints must have getConstraintValue(xarray) getJacobianParticles(xarray) getDotJacobianParticles(varray) methods

/*
All forces has the signature method getForceArray(xarray, varray, masses)
computeExternalForces has signature computeExternalForces(forces,xarray,varray,masses)
//where forces is an array of force objects, xarray is the position array, varray is the velocity array and masses is the masses array
*/
const DEBUG = true;

if (DEBUG){
    console.debug("solver.js loaded","PENDING MODIFY HANDLING OF CONTACT CONSTRAINTS, NEEDED AN ADDITIONAL ITERATION AFTER CLIPPING LAGRANGE MULTIPLIERS (ACTIVE SET METHOD)");
}



//lets write a State class, that will wrap the basic state object and will have some utility methods
//lets make all methods have explicit arguments, so we can use them in the main loop

class StateUtils{
    constructor(){

    }

    unpackState(STATE){
        //lets call the unpackState function
        return stepUnpackState(STATE);
    }
    getForcesObjs(STATE){
        let [_1, forceObjs, _2] = this.unpackState(STATE);
        return forceObjs;
    }
    getConstraintsObjs(STATE){
        let [constraintObjs, _1, _2] = this.unpackState(STATE);
        return constraintObjs;
    }

    getNormalizedConstraintViolations(STATE){
        //we get the constraintObjs and for the distance pin constraints we get the normalized violation by calling the getNormalizedConstraintViolation method
        //that gets as argument the xarray

        let cObjs = this.getConstraintsObjs(STATE);
        let cObjsDistPin = cObjs.filter((c) => c.getType() === "distance" || c.getType() === "pin");
        let normErrors =cObjsDistPin.map((c) => c.getNormalizedConstraintViolation(STATE.xs));

        return normErrors;

    }
    getPolygonsObjs(STATE){
        let [_1, _2, polygonsObjs] = this.unpackState(STATE);
        return polygonsObjs;
    }
    getJacobianMatrices(STATE){
        let s = STATE;
        let [J,dotJ] = computeJacobians(this.getConstraintsObjs(s), s.xs, s.vs);
        return [J,dotJ];
    }
    getExternalForces(STATE){
        let s = STATE;
        let forces = computeExternalForces(this.getForcesObjs(s), s.xs, s.vs, s.masses);
        return forces;
    }
    solveConstraints(STATE, alpha = 0.00, beta = 0.00, clipNegativeContact = false){
        let s = STATE;
        let [J,dotJ] = this.getJacobianMatrices(s);
        let forces = this.getExternalForces(s);
        let out_dict = solveConstraints(J, dotJ, s.xs, s.vs, forces, s.masses, this.getConstraintsObjs(s), alpha, beta, clipNegativeContact);
        return out_dict;
    }

    getAngleBtwParticlesHorizontal(STATE, idxP1, idxP2, unit = "rad"){
        //angle between r_{12} = r_{2} - r_{1} and the horizontal axis [1,0]
        let xs = STATE.xs;
        let [x1,y1] = xs[idxP1];
        let [x2,y2] = xs[idxP2];
        let dr = [x2-x1,y2-y1];
        return getAngle([dr],[1,0],unit)[0];
    }
    getAngleBtwParticlesPairs(STATE,pair1,pair2){
        //pair1 = [idxP1,idxP2], pair2 = [idxP3,idxP4] we'll get r_{12} and r_{34} and get the angle between them
        let xs = STATE.xs;
        let [idxP1,idxP2] = pair1;
        let [idxP3,idxP4] = pair2;
        let [x1,y1] = xs[idxP1];
        let [x2,y2] = xs[idxP2];
        let [x3,y3] = xs[idxP3];
        let [x4,y4] = xs[idxP4];
        let dr12 = [x2-x1,y2-y1];
        let dr34 = [x4-x3,y4-y3];
        return getAngle([dr12],[dr34])[0];
    }

    getAngleBtwConstraints(STATE, idxC1, idxC2, unit = "rad"){
        //constraints_distance is a [[idx1,idx2,distance],...] 
        let constraints = this.getConstraintsObjs(STATE);
        let cs = [constraints[idxC1],constraints[idxC2]];

        // so if the constraint is distance we get idx1 and idx2 and calculate r_{12}
        // if the constraint is pin we get idx and pin_point and calculate r_{idx,pin_point}
        //we then calculate the angle between the vectors
        let xs = STATE.xs;

        let anglesHorizontal = [];//we then sustract angles[1]-angles[0] to get the angle between the constraints
        for (let c of cs){
            if (c.getType() === "distance"){
                let idx1 = c.getParticleIndices()[0];
                let idx2 = c.getParticleIndices()[1];
                let angle = this.getAngleBtwParticlesHorizontal(STATE,idx1,idx2,unit);
                anglesHorizontal.push(angle);
            }
            else if (c.getType() === "pin"){
                let idx = c.getParticleIndices()[0];
                let point = xs[idx];
                let pin_point = c[1]
                let angle = getAngle([math.subtract(point,pin_point)],[1,0],unit)[0];
                anglesHorizontal.push(angle);
            }
            else{
                throw new Error(`Constraint type ${c.getType()} not pin nor distance, but ${c.getType()}`);
            }
        }
        return anglesHorizontal[1]-anglesHorizontal[0];


    }


    getAngleBtwConstraintsDistanceHorizontal(STATE, unit = "rad"){
        //constraints_distance is a [[idx1,idx2,distance],...] array, so we loop over it and get 
        //the angle between r_{12} and the horizontal axis
        let constraints = this.getConstraintsObjs(STATE);
        let xs = STATE.xs;
        let angles = [];
        for (let constraint of constraints){
            if (constraint.getType() === "distance"){
                let idx1 = constraint.getParticleIndices()[0];
                let idx2 = constraint.getParticleIndices()[1];
                let angle = this.getAngleBtwParticlesHorizontal(STATE,idx1,idx2,unit);
                angles.push(angle);
            }
        }
         
        return angles;
    }

    getAngleBtwConstraintsPinHorizontal(STATE,unit = "rad"){
        //constraints_pin is a [[idx,pin_point,pin_distance],...] array, so we loop over it and get 
        //the angle between r_{idx,pin_point} and the horizontal axis
        // returns a [nConstraintsPin] array with the angles
        let constraints = this.getConstraintsObjs(STATE);
        let xs = STATE.xs;
        let angles = [];
        for (let constraint of constraints){
            if (constraint.getType() === "pin"){
                let idx = constraint.getParticleIndices()[0];
                let point = xs[idx];
                let pin_point = constraint.getPinPoint();
                let angle = getAngle([math.subtract(point,pin_point)],[1,0],unit)[0];
                angles.push(angle);
            }
        }
         
        return angles;
    }

    getAngleSegmentsHorizontal(STATE,unit ="rad"){
        //basically any pin or distance constraint constitutes a segment
        // lets refer to the segments with the same index as the constraint
        // constraints are [distance,pin,contact]
        // so we'll have the same segments as [distance,pin]
        //this function will return an array of angles between the segments and the horizontal axis
        let anglesDistance = this.getAngleBtwConstraintsDistanceHorizontal(STATE,unit);
        let anglesPin = this.getAngleBtwConstraintsPinHorizontal(STATE,unit);
        return anglesDistance.concat(anglesPin);
    }
    getDotAngleSegments(STATE,unit ="rad"){
        // here we'll return the omega for each segment
        // we dont refer to them as Horizontal because the omega is a time derivative that
        // doesnt depend on the reference theta0
        // so basically for a segment omega = \vec{v__{rel12}}\vec{\theta}/R_{12} where R_{12} is the distance between the particles
        // and \vec{v__{rel12}} is the relative velocity between the particles
        //and \vec{\theta} is the unit vector perpendicular to the segment
        // so we'll return an array of nsegments with the omega for each segment
        //\vec{\theta} = [-dry,drx]/norm([-dry,drx]) = [-dr12[1],dr12[0]]/norm(dr12)

        //lets get the [distance,pin] constraints from the state constraints list
        //we filter by getType (distance,pin)
        //we then get the particle indices and calculate the relative velocity
        // we then calculate the omega for each segment
        let constraints = this.getConstraintsObjs(STATE).filter((c) => c.getType() === "distance" || c.getType() === "pin");
        let omegas = [];

        for (let c of constraints){
            if (c.getType() === "distance"){
                let idx1 = c.getParticleIndices()[0];
                let idx2 = c.getParticleIndices()[1];
                let dr = math.subtract(STATE.xs[idx2],STATE.xs[idx1]);
                let dr_norm = math.norm(dr);
                let vrel = math.subtract(STATE.vs[idx2],STATE.vs[idx1]);
                let omega = math.dot(vrel,[-dr[1],dr[0]])/dr_norm;
                omegas.push(omega);
            }
            else if (c.getType() === "pin"){
                let idx = c.getParticleIndices()[0];
                let dr = math.subtract(STATE.xs[idx],c.getPinPoint());
                let dr_norm = math.norm(dr);
                let vrel = STATE.vs[idx];
                let omega = math.dot(vrel,[-dr[1],dr[0]])/dr_norm;
                omegas.push(omega);
            }
            else{
                throw new Error(`Constraint type ${c.getType()} not pin nor distance, but ${c.getType()}`);
            }

        }

        if (unit === "deg"){
            omegas = omegas.map((omega) => math.unit(omega,"rad").toNumber("deg"));
        }



        return omegas;

        }
   

    getKineticEnergy(STATE){
        let masses = STATE.masses;
        let vs = STATE.vs;
        return calculateKineticEnergy(vs,masses);
    }

    getPotentialEnergy(STATE){
        let g = STATE.gravity[0];
        //we asume direction is [0,-1] so mg*xs[1] is the potential energy for each particle
        let xs = STATE.xs;
        let masses = STATE.masses;
        let n = xs.length;
        let pe = 0;
        for (let i = 0;i<n;i++){
            pe += masses[i]*g*xs[i][1];
        }
        return pe;
    }

    getTotalEnergy(STATE){
        //neglecting the spring potential energy
        return this.getKineticEnergy(STATE)+this.getPotentialEnergy(STATE);
    }

    getCOM(STATE){
        let masses = STATE.masses;
        let xs = STATE.xs;
        return calculateCOM(xs,masses);
    }

    getNeighborsDict(STATE){
        /*
        Returns two dictionaries, neighborsDict and pinDict
        neighborsDict example {0:[[1,distance01],[2,distance02]],1:[[0,distance01],[2,distance12]],2:[[0,distance02],[1,distance12]]}
        pinDict example {0:[pin_point0,pin_distance0],1:[pin_point1,pin_distance1],2:[pin_point2,pin_distance2]}

        */
        let [neighborsDict,pinDict] = getNeighborsDict(STATE);
        return [neighborsDict, pinDict];
    }


}
    




/*
GENERAL FUNCTIONALITIES
*/




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
                            console.debug(`NEW COLLISION :Particle ${i} with position ${point} closest edge ${closest_edge} distance ${distance} threshold ${threshold}`);

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
                        console.debug(`RESOLVED COLLISION:Particle ${particle_index} with position ${xs[particle_index]} edge ${edge_index} distance ${current_distance} treehold ${treshold}`);
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



            let [_1,_2,polygonObjs] = stepUnpackState(STATE);
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
     alpha = 0.00, beta = 0.00, clipNegativeContact = true){


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
    //lets add for stability a small value to the diagonal of A if nconstraints > 0
    
    if (nconstraints > 0){
        let eps = 1e-6;
        let I = math.identity(nconstraints);
        A = math.add(A, math.multiply(eps,I));
    }


    let lagrange_multipliers = math.lusolve(A,b).reshape([-1]);
    // We have to clip the contact lagrange multipliers to be >= 0
    // lets get the indices of the contact constraints using the getType method

    if (clipNegativeContact){
        let contact_indices = constraints.map((constraint,i) => [constraint.getType(),i]).filter(c => c[0] === "contact").map(c => c[1]);
        for (let i = 0; i < contact_indices.length; i++){
            let index = contact_indices[i];
            if (lagrange_multipliers.get([index]) < 0){
                lagrange_multipliers.set([index],0);
            }
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

function stepUnpackState(STATE){
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
        forceObjs = forceObjs.concat(s.damping.map(f => new Damping(f[0], f[1])));
    }
    if (s.hasOwnProperty("springs")){
        forceObjs = forceObjs.concat(s.springs.map(f => new Spring(f[0], f[1], f[2], f[3])));
    }
    if (s.hasOwnProperty("springs2points")){
        forceObjs = forceObjs.concat(s.springs2points.map(f => new Spring2Point(f[0], f[1], f[2], f[3])));
    }




    return [constraintsObjs, forceObjs, polygonsObjs];

}
    



function stepSemiEulerActiveSet(STATE,PARAMETERS, callbacksActuators = null){
    //There is an issue that becomes apparent when particles in contact with a polygon detach from it, I think it stems in part from the handling of the inequality constraints
    // in the solver so far, we handle the inequality constraints by clipping the lagrange multipliers to be >= 0, the thing is, even if we prevent contact forces to be attractive
    // the rest of the constraint forces are calculated as if the contact force were a pin constraint, so we need to to do a following interation after clipping the lagrange multipliers
    // we need to remove the contact constraints with negative lagrange multipliers and solve the system again, repeating the process if a contact lagrange multiplier becomes negative after removing a contact constraint.
    // This is a specific instance of the active set method.

    //I'll try to describe all of this in more detail

    /*
    * We have to find the lagrange multipliers, get the most negative contact lagrange multipliers 
    and we remove it from the active set we solve now the system again but with tne updated active set
    What we'll do is remove that constraint and recalculate jacobians etc and solve again the reduced constraints system

    */
    
    /*
    * 1 Unpack the state object
    * 2 Calculate external forces
    * 3 Calculate constraint forces
    *   If contact lag multipliers are<0 their respective collision and constraintObj is removed
    *  We solve the system again and repeat the process
    * Now, there is one thing, one removing a given contact constraint, the previous removed constraint can become active again
    * This might be specially problematic for equilibrium states with multiple contacts for bodies that are mutually equilibrated
    * By now we'll just remove the constraint and solve the system again, in the future we must do line search and check wether
    * some normal negative acceleration appear for some of the constraints that were previously removed

    */

    //callbackActuator is a temporary solution
    // to implement feedback control, all the callbacks will be called with the state object
    // actuatorForce = callbackActuator(state) -> [nparticles,2]
    // we just add this force to the external forces
    let mandatoryParameters = ["collisionThreshold", "muFriction", "alpha", "beta"];
    for (let key of mandatoryParameters){
        if (!PARAMETERS.hasOwnProperty(key)){
            throw new Error(`The PARAMETERS object must have the key ${key}`)
        }
    }
    let s = STATE;

    //lets save the original collisions
    let originalCollisions = JSON.parse(JSON.stringify(STATE.collisions));
    let [constraintObjs, forceObjs, polygonsObjs] = stepUnpackState(STATE);
    let externalForces = computeExternalForces(forceObjs, s.xs, s.vs, s.masses);
    let actuatorForces = math.zeros(s.xs.length,2);
    if (callbacksActuators !== null){
        for (let callbackActuator of callbacksActuators){
            let actuatorForce = callbackActuator(STATE);
            actuatorForces = math.add(actuatorForces, actuatorForce);
            
        }

    }
    externalForces = math.add(externalForces, actuatorForces);
    //collisions are an array of [iP,jE,Pidx] where iP is the particle index, jE is the edge index and Pidx is the polygon index

    function solveConstraintsAndFindInactive( STATE, PARAMETERS, externalForces){
        let [constraintObjs, forceObjs, polygonsObjs] = stepUnpackState(STATE);
    
        let s = STATE;
        //collisions are an array of [iP,jE,Pidx] where iP is the particle index, jE is the edge index and Pidx is the polygon index
        let [J,dotJ] = computeJacobians(constraintObjs, s.xs, s.vs);
        let outSolverConstraints = solveConstraints(J, dotJ, s.xs, s.vs, externalForces,
             s.masses, constraintObjs, PARAMETERS.alpha, PARAMETERS.beta, false);
        let indicesConstraintsContact = constraintObjs.map((c,i) => [c.getType(),i]).filter(c => c[0] === "contact").map(c => c[1]);
        let lagrangeMultipliersContact = indicesConstraintsContact.map(i => outSolverConstraints.lagrange_multipliers[i]);
        let indicesNegative = lagrangeMultipliersContact.map((l,i) => [l,i]).filter(l => l[0] < 0).map(l => l[1]);
        //lets print the lagrange multipliers, and the indicesNegative
        //console.debug("Lagrange multipliers contact", lagrangeMultipliersContact);
        //console.debug("Indices negative", indicesNegative);

        return [outSolverConstraints, indicesNegative, J, dotJ];
    }

    let whileCondition = true;
    let [outSolverConstraints, indicesNegative,J,dotJ] = [null,[],null,null];
    while (whileCondition){
        [outSolverConstraints, indicesNegative,J,dotJ] = solveConstraintsAndFindInactive(STATE, PARAMETERS, externalForces);
        if (indicesNegative.length === 0){
            whileCondition = false;
        }else{
            // we just update the state by removing the collision with the most negative lagrange multiplier
            let activeCollisions = s.collisions.filter((c,i) => !indicesNegative.includes(i));
            //console.debug(`Removing collision ${indicesNegative[0]} with lagrange multiplier ${outSolverConstraints.lagrange_multipliers[indicesNegative[0]]}`);
            STATE.collisions = activeCollisions;
        }
    }

    [constraintObjs, forceObjs, polygonsObjs] = stepUnpackState(STATE);
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

    //lets clip the velocities for stability to < 6 in module
    newVs = newVs.map(v => {
        let mod = math.norm(v);
        if (mod > 4){
            return math.multiply(4/mod,v);
        }
        return v;
    }); 
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
        "actuator_forces": actuatorForces.toArray(),
        "total_forces": totalForces,
        "lagrange_multipliers": outSolverConstraints.lagrange_multipliers,
        "lagrange_multipliers_contact": lagrangeMultipliersContact,
        "constraints_vals": outSolverConstraints.constraints_vals,
        "collisions": STATE.collisions,
        "J": J.toArray(),
        "dot_J": dotJ.toArray(),

    }
    let [collisionsArray, newCollisionsArray, resolvedCollisionsArray] = colHandler.updateCollisions(STATE, PARAMETERS.collisionThreshold);
    
    newState.collisions = collisionsArray;
    newState.new_collisions = newCollisionsArray;
    newState.resolved_collisions = resolvedCollisionsArray;
    //each collision entry is [iP,jE,Pidx] where iP is the particle index, jE is the edge index and Pidx is the polygon index

    //lets add to newState the state properties not included in the previous object
    //we loop over STATE keys and if they are not present in newState we add them
    for (let key in STATE){
        if (!newState.hasOwnProperty(key)){
            newState[key] = STATE[key];
        }
    }

    return newState;
}



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
    let [J, dot_J] = computeJacobians(constraintsObjs, s.xs, s.vs);
    let outSolveConstraints = solveConstraints(J, dot_J, s.xs, s.vs,
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
    s["constraints_vals"] = outSolveConstraints.constraints_vals;
    s["lagrange_multipliers"] = outSolveConstraints.lagrange_multipliers;
    s["J"] = J.toArray();
    s["dot_J"] = dot_J.toArray();
    s["time"] = 0;

    return s;



    }


function getConstraintsRigid(xs){
    //The distance constraints that make the system a rigid body, from delauney triangulation
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

    return constraints_distance;

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





function _postProcessStateStory(stateStory){
    //lets add the kinetic energy and the angle between the pin and the horizontal
    let sUtils = new StateUtils();
    let nsteps = stateStory.length;
    for(let i = 0; i < nsteps; i++){
        let state = stateStory[i];
        state.kineticEnergy = sUtils.getKineticEnergy(state);
        state.totalEnergy = sUtils.getTotalEnergy(state);
        state.COM = sUtils.getCOM(state);
        state.angleConstraintsPin = sUtils.getAngleBtwConstraintsPinHorizontal(state, "deg")[0];
        state.angleConstraintsDistance = sUtils.getAngleBtwConstraintsDistanceHorizontal(state, "deg")[0];
    }
}   

function integrateSystem(STATE, nsteps, callbacksActuators = null,
                     method = "SemiEulerActiveSet",
                        PARAMETERS = {})
                        {
    

    let defaultParams = {
        "alpha": 45,
        "beta": 45,
        "dt": 0.001,
        "muFriction": 0,
        "collisionThreshold": 0.0001,
    }

    for (let key in defaultParams){
        if (!PARAMETERS.hasOwnProperty(key)){
            console.debug(`PARAMETERS object does not have key ${key}, using default value ${defaultParams[key]}`);
            PARAMETERS[key] = defaultParams[key];
        }
    }

    let s = initState(STATE, PARAMETERS);

    let STATE_STORY = [s];
    for (let i = 0; i < nsteps; i++){
        //lets copy the state object, just in case
        let state = JSON.parse(JSON.stringify(s));
        if (method === "SemiEuler"){
            let newState = stepSemiEuler(state,PARAMETERS);
            STATE_STORY.push(newState);
            s = newState;
        }//else if method  SemiEulerActiveSet
        
        else if (method === "SemiEulerActiveSet"){
            let newState = stepSemiEulerActiveSet(state,PARAMETERS, callbacksActuators);
            STATE_STORY.push(newState);
            s = newState;
             }

        else {
            throw new Error(`Method ${method} not implemented`)
        }
    }

    _postProcessStateStory(STATE_STORY);
    return STATE_STORY;
    }



 function  _checkInitState(STATE){
        //here we'll write zeros for the missing properties in the state
        // and print a info message to the console
        // we'll check for the forces, if not we init them to [nparticles,2] array of zeros
        //we'll check for polygons, if not we init them to an empty array
        //We'll check for constraints and if not we init them to an empty array
        //vs to zeros if not present aswell
        //xs and masses are mandatory

        if (!STATE.hasOwnProperty("xs")){
          throw new Error("Missing xs in STATE");
        }
        if (!STATE.hasOwnProperty("masses")){
          throw new Error("Missing masses in STATE");
        }

        
        let nparticles = STATE.xs.length;
        let names2InitZerosIfMissing = ["vs","external_forces", "constraint_forces",
                                           "contact_forces", "friction_forces",
                                           "actuator_forces",
                                            "total_forces"];
        let names2InitEmptyArrayIfMissing = ["polygons", "constraints_distance",
                                             "constraints_pin","collisions","new_collisions","resolved_collisions",
                                             "springs","springs2points"];

        for (let name of names2InitZerosIfMissing){
          if (!STATE.hasOwnProperty(name)){
            STATE[name] = Array(nparticles).fill().map(()=>Array(2).fill(0));
            console.info(`Property ${name} not found in STATE, initializing to zeros`);
          }
        }

        for (let name of names2InitEmptyArrayIfMissing){
          if (!STATE.hasOwnProperty(name)){
            STATE[name] = [];
            console.info(`Property ${name} not found in STATE, initializing to empty array`);
          }
        }


        if (!STATE.hasOwnProperty("time")){
          STATE["time"] = 0;
          console.info("Property time not found in STATE, initializing to 0");
        }


        return STATE;
      }

function step(STATE, PARAMETERS = {}, callbacksActuators = null){

    //updates the state object with the new state after a time step

    _checkInitState(STATE);


    let defaultParams = {
        "alpha": 35,
        "beta":35,
        "dt": 0.01,
        "muFriction": 0.25,
        "collisionThreshold": 0.0001,

    }

    for (let key in defaultParams){
        if (!PARAMETERS.hasOwnProperty(key)){
            console.debug(`PARAMETERS object does not have key ${key}, using default value ${defaultParams[key]}`);
            PARAMETERS[key] = defaultParams[key];
        }
    }

    let newState = stepSemiEulerActiveSet(STATE,PARAMETERS, callbacksActuators);

    //lets loop over newState keys and update the STATE object, so the outer scope object is updated
    for (let key in newState){
        STATE[key] = newState[key];
    }


}




export { stepSemiEulerActiveSet, calculateCOM, calculateKineticEnergy,
     integrateSystem,initState, calculateInertiaMoment,
        computeJacobians, solveConstraints, stepUnpackState,
         collisions2ContactConstraints,computeExternalForces, calculateConstraintForces, calculateConstraintForcesFromTrajectory,
         CollisionHandler, ConstraintContact, Polygon, StateUtils, getConstraintsRigid,
         step
        };