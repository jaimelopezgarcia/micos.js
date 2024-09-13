class ConstraintDistance {
    constructor(particle1_idx, particle2_idx, distance) {
      this.idx1 = particle1_idx;
      this.idx2 = particle2_idx;
      this.distance = distance;
    }

    getType() {
      return "distance";
    }
  
    getParticleIndices() {
      return [this.idx1, this.idx2];
    }
    
    
    distfunc(x1, y1, x2, y2) {
        return Math.sqrt((x2-x1)**2 + (y2-y1)**2);
    }
    getDistanceParticles(xarray) {
    // lets throw an error if xarray 2 dimension is not 2

    if (xarray[0].length != 2){
        throw new Error("xarray must be a 2D array");
    }

      // Calculate the distance between the two particles
        let x1 = xarray[this.idx1][0];
        let y1 = xarray[this.idx1][1];
        let x2 = xarray[this.idx2][0];
        let y2 = xarray[this.idx2][1];
        return this.distfunc(x1, y1, x2, y2);
    }


  
    getConstraintValue(xarray,squared = true) {
      // C(x) = current_distance(x)**2-distance**2
      if(squared){
        return this.getDistanceParticles(xarray)**2 - this.distance**2;
      }
      else{
        return this.getDistanceParticles(xarray) - this.distance;
      }
    }

    getNormalizedConstraintViolation(xarray){
      //we give the non squared constraint value divided by the distance
      let val = this.getConstraintValue(xarray, false);
      let dist = this.getDistanceParticles(xarray);
      return val/dist;
    } 
  
    getJacobianParticles(xarray){
      // Jac of the distance constraint for  the two particles, it must be a (2,2) array 
      // with the gradient of the distance constraint for each particle
        let x1 = xarray[this.idx1][0];
        let y1 = xarray[this.idx1][1];
        let x2 = xarray[this.idx2][0];
        let y2 = xarray[this.idx2][1];

        let J = [[x1-x2, y1-y2],[x2-x1, y2-y1]]
  
      return J;
    }
  
    getDotJacobianParticles(varray){
      // Jac of the time derivative of the distance constraint for  the two particles, it must be a (2,2) array
      // with the time derivative of J, basically the v^{t}Hessian_constraint
        let v1x = varray[this.idx1][0];
        let v2x = varray[this.idx2][0];
        let v1y = varray[this.idx1][1];
        let v2y = varray[this.idx2][1];
        let Jdot = [[v1x-v2x, v1y-v2y],[v2x-v1x, v2y-v1y]]
  
      return Jdot;
    }
}

class ConstraintPin {
    constructor(particle_idx, pin_point, distance) {
        // lets throw an error if pin_point is not a 2D array
        if (pin_point.length != 2){
            throw new Error("pin_point must be a 2D array");
        }
      this.idx = particle_idx;
      this.pin_point = pin_point;
      this.distance = distance;
    }

    getType() {
        return "pin";
    }
    
    getPinPoint() {
        return this.pin_point;
    }
    getParticleIndices() {
      return [this.idx];
    }
    
    distfunc(x1, y1, x2, y2) {
        return Math.sqrt((x2-x1)**2 + (y2-y1)**2);
    }


    getDistancePoint(xarray) {
    // lets throw an error if xarray 2 dimension is not 2
    if (xarray[0].length != 2){
        throw new Error("xarray must be a 2D array");
    }

      // Calculate the distance between the two particles
        let x1 = xarray[this.idx][0];
        let y1 = xarray[this.idx][1];
        let x2 = this.pin_point[0];
        let y2 = this.pin_point[1];
        return this.distfunc(x1, y1, x2, y2);
    }
    getConstraintValue(xarray, squared = true) {
        // C(x) = current_distance(x)**2-distance**2
        if (squared){
            return this.getDistancePoint(xarray)**2 - this.distance**2;
        }
        else{
            return this.getDistancePoint(xarray) - this.distance;
        }
        }

    getNormalizedConstraintViolation(xarray){
      //we give the non squared constraint value divided by the distance
      let val = this.getConstraintValue(xarray, false);
      let dist = this.getDistancePoint(xarray);
      return val/dist;
    }

    getJacobianParticles(xarray){
        // Jac of the distance constraint for  the two particles, it must be a (2,2) array 
        // with the gradient of the distance constraint for each particle
        let x1 = xarray[this.idx][0];
        let y1 = xarray[this.idx][1];
        let x2 = this.pin_point[0];
        let y2 = this.pin_point[1];
        let J = [[x1-x2, y1-y2]]
        return J;
    }

    getDotJacobianParticles(varray){
        // Jac of the time derivative of the distance constraint for  the two particles, it must be a (2,2) array
        // with the time derivative of J, basically the v^{t}Hessian_constraint
        let v1x = varray[this.idx][0];
        let v1y = varray[this.idx][1];
        let Jdot = [[v1x, v1y]]
        return Jdot;
    }
}


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

export { ConstraintDistance, ConstraintPin, ConstraintContact};