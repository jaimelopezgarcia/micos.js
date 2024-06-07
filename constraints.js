class ConstraintDistance {
    constructor(particle1_idx, particle2_idx, distance) {
      this.idx1 = particle1_idx;
      this.idx2 = particle2_idx;
      this.distance = distance;
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


  
    getConstraintValue(xarray) {
      // C(x) = current_distance(x)-distance
      return this.getDistanceParticles(xarray) - this.distance;
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
    getConstraintValue(xarray) {
        // C(x) = current_distance(x)-distance
        return this.getDistancePoint(xarray) - this.distance;
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

export { ConstraintDistance, ConstraintPin };