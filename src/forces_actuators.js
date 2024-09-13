import {getNeighborsDict,getAngle} from "./math_utils.js";
//lets import math.js
import * as math from "mathjs";

/*
All forces has the signature method getForceArray(xarray, varray, masses)
computeExternalForces has signature computeExternalForces(forces,xarray,varray,masses)
//where forces is an array of force objects, xarray is the position array, varray is the velocity array and masses is the masses array
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
























/////////
/*
ACTUATORS
*/
/////////


class ForceActuator{
    //just a external force that can be applied to selected particles
    constructor(particle_ids, force){
        this.particle_ids = particle_ids;
        this.force = force;//[2] array
        
    }
    getName(){
        return "ForceActuator";
    }
    getForceArray(STATE){

        // forceArray must be a [nparticles,2] array
        //so we initialize it with zeros
        let forceArray = math.zeros([STATE.xs.length,2]);

        for(let i = 0; i < this.particle_ids.length; i++){
            let particle_id = this.particle_ids[i];
            forceArray[particle_id] = this.force;
            
            

        }
        console.log
        return forceArray
        

    }
}









//lets create a linear actuator, basically a spring that changes its rest length


class SpringActuator{
    /*
    Example usage SpringActuator:
    let contractionFunction = (distance, time, l0Current, kCurrent) => [Math.cos(time)+1, kCurrent];
    let actuator = new SpringActuator(0, 1, contractionFunction, 1, 0.5);
    let callbacksActuators = [actuator.getForceArray.bind(actuator)];
    */
    constructor(idxParticle1, idxParticle2,
        contractionFunction, springConstant = 1, restLength = 0.5){
        //contractionFunction is a function that takes as argument distance betwen parts 1 and 2, time,current l0 and k and returns the new rest length and spring constant
        this.idx1 = idxParticle1;
        this.idx2 = idxParticle2;
        this.k = springConstant;
        this.l0 = restLength;
        this.spring = new Spring(idxParticle1, idxParticle2, springConstant, restLength);//we change the restLength later
        this.contractionFunction = contractionFunction;
        this._checkContractionFunction();
    }

    _checkContractionFunction(){
        // basic check to see the contractionFunction has the right signature and gives a number
        let [distance, time,l0Current, kCurrent] = [1,1,1,1];
        let [l0, k] = this.contractionFunction(distance, time, l0Current, kCurrent);
        if(typeof l0 != "number" || typeof k != "number"){
            throw "The contraction function must return a number";
        }
    }

    getName(){
        return `SpringActuator_${this.idx1}_${this.idx2}`;
    }


    getForceArray(STATE){

        let [xarray,varray,time] = [STATE.xs, STATE.vs, STATE.time];
        let x1 = xarray[this.idx1];
        let x2 = xarray[this.idx2];
        let x1_x = x1[0];
        let x1_y = x1[1];
        let x2_x = x2[0];
        let x2_y = x2[1];
        let distance = Math.sqrt((x2_x-x1_x)**2 + (x2_y-x1_y)**2);
        let [l0, k] = this.contractionFunction(distance, time, this.l0, this.k);
        //lets throw an error if l0 or k are not numbers or they are negative
        if(typeof l0 != "number" || typeof k != "number" || l0 < 0 || k < 0){
            throw "The contraction function must return positive numbers";
        }
        this.l0 = l0;
        this.k = k;
        this.spring.rest_length = l0;
        this.spring.spring_constant = k;
        return this.spring.getForceArray(xarray, varray, STATE.masses);

    }

    getForce(STATE){
        //lets return the a [2,2] array with the forces acting on the particles
        let forceArray = this.getForceArray(STATE);
        let f1 = forceArray[this.idx1];
        let f2 = forceArray[this.idx2];
        return [f1, f2];

}


}



function calculateTorqueLeastNormedForce(pivotPoint, applicationPoint, torqueValue){
    // here we calculate the least normed force that can generate the given torque (perpendicular to the r vector)
    // because it is a 2D problem we can use for the torque vector [0,0,torqueValue]
    //and we pick the force vector [fx,fy,0]

    //returns the force vector [fx,fy] that generates the torqueValue at the
    // applicationPoint respect to the pivotPoint

    let torqueVec = [0,0,torqueValue];
    let rVec = math.subtract(applicationPoint, pivotPoint);
    let normR = math.norm(rVec);
    let rVecUnit = math.divide(rVec, normR);
    //lets make rvec a 3D vector
    rVecUnit.push(0);
    let forceVec = math.divide( math.cross(torqueVec,rVecUnit),normR);
    //we retrieve the first 2 components of the force vector
    return forceVec.slice(0,2);
}



class TorqueActuatorPin{
    /*
Example usage:
let torqueFun = (theta, omega, time) => torqueFunction(theta, omega, time, 0.0, 50, 0.5);
let actuator = new TorqueActuatorPin([0,0], 0, torqueFun); pinpoint targetParticletorqueFunction
let callbacksActuators = [actuator.getForceArray.bind(actuator)];

    */
        constructor(pinPoint, targetParticle, torqueFunction){
            this.pinPoint = pinPoint;
            this.targetParticle = targetParticle;
            this.torqueFunction = torqueFunction;
            this._checkTorqueFunction(torqueFunction);
        }

        getName(){
            //we return TorqueActuatorPin_targetParticleIdx
            return `TorqueActuatorPin_${this.targetParticle}`;
        }
        _checkTorqueFunction(torqueFunction){
            // basic check to see the torqueFunction has the right signature and gives a number
            let [theta, omega, time] = [0,0,0];
            let torque = torqueFunction(theta, omega, time);
            if(typeof torque != "number"){
                throw "The torque function must return a number";
            }
        }

        getForceArray(STATE, returnTorque = false){
            let [xarray,varray,time] = [STATE.xs, STATE.vs, STATE.time];
            let targetIdx = this.targetParticle;
            let rpivotTarget = math.subtract(xarray[targetIdx], this.pinPoint);
            let thetaTarget = getAngle(rpivotTarget, [1,0]);
            let [rTemp, vTemp] = [[rpivotTarget[0], rpivotTarget[1],0], [varray[targetIdx][0], varray[targetIdx][1],0]];
            let omegaTarget = math.cross(rTemp, vTemp)[2]/math.norm(rpivotTarget);
            let torque = this.torqueFunction(thetaTarget, omegaTarget, time);
            let forceTarget = calculateTorqueLeastNormedForce(this.pinPoint, xarray[targetIdx], torque);
            let forceArray = Array(xarray.length).fill().map((_,i) => [0,0]);
            forceArray[targetIdx] = forceTarget;
            
            if(returnTorque){
                return [forceArray, torque];
            }else{
                return forceArray;
            }
       
        }

    getTorque(STATE){
        //util function to retrieve the torque associated with the pinpoint target joint
        // basically we call getForceArray and return the second element
        return this.getForceArray(STATE, true)[1];
    }

    }

class TorqueActuatorJoint{
    constructor(pivotParticleIdx, targetParticleIdx, torqueFunction){
        this.pivotParticleIdx = pivotParticleIdx;
        this.targetParticleIdx = targetParticleIdx;
        this.torqueFunction = torqueFunction;
        this.firstCall = true;// on the first getForceArray call we assert that pivot and target particles are neighbors
        
        this._checkTorqueFunction(torqueFunction);
    }

    getName(){
        //we return TorqueActuatorJoint_pivotParticleIdx_targetParticleIdx
        return `TorqueActuatorJoint_${this.pivotParticleIdx}_${this.targetParticleIdx}`;
    }
    _checkTorqueFunction(torqueFunction){
        // basic check to see the torqueFunction has the right signature and gives a number
        let [theta, omega, time] = [0,0,0];
        let torque = torqueFunction(theta, omega, time);
        if(typeof torque != "number"){
            throw "The torque function must return a number";
        }
    }

    _checkNeighbors(STATE){
        if (this.firstCall){
                let [neighborsDict, pinDict] = getNeighborsDict(STATE);
                this.neighborsDict = neighborsDict;
                this.pinDict = pinDict;
                this.firstCall = false;

                //lets check that the pivot and target particles are neighbors
                let neighborsPivot = this.neighborsDict[this.pivotParticleIdx];
                //lets throw an error if target particles is not in neighborsPivot
                let targetParticleInNeighbors = neighborsPivot.some(([idx, _]) => idx == this.targetParticleIdx);
                if(!targetParticleInNeighbors){
                    throw `The target particle ${this.targetParticleIdx} is not a neighbor of the pivot particle ${this.pivotParticleIdx}, neighbors are ${neighborsPivot}`;
                }
    
            }
    }
    getForceArray(STATE, returnTorque = false){

        this._checkNeighbors(STATE);

        let [xarray,varray,time] = [STATE.xs, STATE.vs, STATE.time];
        let pivotXs = xarray[this.pivotParticleIdx];

        let torqueActuatorPin = new TorqueActuatorPin(pivotXs, this.targetParticleIdx, this.torqueFunction);

        return torqueActuatorPin.getForceArray(STATE, returnTorque);
    }

    getTorque(STATE){

        return this.getForceArray(STATE, true)[1];
    }

}


// lets make an elbow joint actuator to implement a more physical torque actuator
//this will only apply on particles/joints that have atleast 2 neighbors
// we'll make sure that net angular momentum  and linear momentum is conserved
// for this we impart opposing torques on the particle neighbors, and we apply an additional force on the joint so F_joint+ F_neighbor1 + F_neighbor2 = 0
//we'll throw an error if the targeted particle to act as a joint has other than 2 neighbors.

class TorqueActuatorElbowJoint{

    constructor(jointParticleIdx, targetParticleIdx, torqueFunction){
        this.jointParticleIdx = jointParticleIdx;
        this.torqueFunction = torqueFunction;
        this.firstCall = true;// on the first getForceArray call we assert that joint particle has 2 neighbors
        //we'll impart a positive torque on the targetParticleIdx and negative torque on the other neighbors
        this.targetParticleIdx = targetParticleIdx;

        
        this._checkTorqueFunction(torqueFunction);
    }

    getName(){
        return `TorqueActuatorElbowJoint_${this.jointParticleIdx}`;
    }

    _checkTorqueFunction(torqueFunction){
        // basic check to see the torqueFunction has the right signature and gives a number
        let [theta, omega, time] = [0,0,0];
        let torque = torqueFunction(theta, omega, time);
        if(typeof torque != "number"){
            throw "The torque function must return a number";
        }
    }

    _checkNeighbors(STATE){
        if (this.firstCall){
                let [neighborsDict, pinDict] = getNeighborsDict(STATE);
                this.neighborsDict = neighborsDict;
                this.pinDict = pinDict;
                this.firstCall = false;

                //lets check that the joint particle has 2 neighbors and the target particle is one of them
                let neighborsJoint = this.neighborsDict[this.jointParticleIdx];
                //lets throw an error if joint particle has other than 2 neighbors
                if(neighborsJoint.length != 2){
                    throw `The joint particle ${this.jointParticleIdx} has other than 2 neighbors, neighbors are ${neighborsJoint}`;
                }
                //lets throw an error if target particles is not in neighborsJoint
                let targetParticleInNeighbors = neighborsJoint.some(([idx, _]) => idx == this.targetParticleIdx);

                if(!targetParticleInNeighbors){
                    throw `The target particle ${this.targetParticleIdx} is not a neighbor of the joint particle ${this.jointParticleIdx}, neighbors are ${neighborsJoint}`;
                }


    
                }
             }

    getForceArray(STATE, returnTorque = false){
        this._checkNeighbors(STATE);
        let [xarray,varray,time] = [STATE.xs, STATE.vs, STATE.time];
        let neighborsJoint = this.neighborsDict[this.jointParticleIdx];
        let [idx1, idx2] = neighborsJoint.map(([idx, _]) => idx);
        let [idxTarget, idxOther] = neighborsJoint.find(([idx, _]) => idx == this.targetParticleIdx) ? [this.targetParticleIdx, idxOther] : [idxOther, this.targetParticleIdx];
        let pivotXs = xarray[this.jointParticleIdx];
        let torqueActuatorTarget = new TorqueActuatorPin(pivotXs, idxTarget, this.torqueFunction);
        //on the other neighbor we apply the opposite torque
        let [forceArrayTarget,torqueTarget] = torqueActuatorTarget.getForceArray(STATE,true);
        //lets create a torqueActuator for the other neighbor with a custom torque function that gives -torqueTarget
        let torqueFunctionOther = (theta, omega, time) => -torqueTarget;
        let torqueActuatorOther = new TorqueActuatorPin(pivotXs, idxOther, torqueFunctionOther);
        let [forceArrayOther,torqueOther] = torqueActuatorOther.getForceArray(STATE,true);
        let forceArrayJoint = math.multiply(-1, math.add(forceArrayTarget, forceArrayOther));
        let forceArray = math.add(forceArrayJoint, math.add(forceArrayTarget, forceArrayOther));
        if(returnTorque){
            return [forceArray, torqueTarget];
        }
        return forceArray;
    }

    getTorque(STATE){
        return this.getForceArray(STATE, true)[1];
    }

}









export {ConstantForce,Gravity,Damping,Spring,Spring2Point,computeExternalForces, calculateTorqueLeastNormedForce,
    TorqueActuatorPin, TorqueActuatorJoint, SpringActuator, ForceActuator, TorqueActuatorElbowJoint
}