/*
Requires math.js to be loaded on the page
*/

function getAngle(vectorsQuery, vectorTarget, angle_unit="rad"){
    //vectorsQuery is a [npoints, 2] array of vectors whose angle with respect to vectorTarget

    //we'll calculate. This function must return a [npoints] array of angles
    
    if (!Array.isArray(vectorsQuery) || !Array.isArray(vectorTarget)){
        throw new Error(`Invalid type of vectorsQuery or vectorTarget. Expected Array, got ${typeof vectorsQuery} and ${typeof vectorTarget}`);
    }

    if ( getArrayDimensions(vectorsQuery) > 2){
        throw new Error(`Invalid dimensions of vectorsQuery. Expected 1 or 2, got ${getArrayDimensions(vectorsQuery)}`);
    }

    // if vectorsQuery is a single vector we'll [vectorQuery]
    if (getArrayDimensions(vectorsQuery) == 1){
        vectorsQuery = [vectorsQuery];
    }
    
    let angles = vectorsQuery.map(v => math.unit(math.atan2(v[1], v[0]) - math.atan2(vectorTarget[1],
                                         vectorTarget[0]), "rad").toNumber(angle_unit));
    
    return angles;
}




function getArrayShape(array) {
    // example usage getArrayShape([ [1,2],[3,4] ]) returns [ 2, 2 ]
    let shape = [];
    while (Array.isArray(array)) {
      shape.push(array.length);
      array = array[0];
    }
    return shape;
  }

  function getArrayDimensions(array) {
    //example usage getArrayDimensions([[[1,2],[3,4]],[[5,6],[7,8]]]) returns 3
    let shape = getArrayShape(array);
    return shape.length;
    }

function isShapeEqual(array1, array2) {
    let shape1 = getArrayShape(array1);
    let shape2 = getArrayShape(array2);
    return JSON.stringify(shape1) === JSON.stringify(shape2);
    }





function translate(points, shift){
    //points is a [npoints, 2] array
    //shift is a [2] array
    //lets throw an error if points or shift are not arrays
    if (!Array.isArray(points) || !Array.isArray(shift)){
        throw new Error(`Invalid type of points or shift. Expected Array, got ${typeof points} and ${typeof shift}`);
    }
    // lets throw an error on dimensions mismatch of points and shift
    if (points[0].length != 2 || shift.length != 2){
        throw new Error("Invalid dimensions of points or shift");
    }
    let p_T = points.map(p => math.add(p, shift));
    return p_T;
}

function rotate(points, angle, center, angle_unit="deg"){
    //points is a [npoints, 2] array
    //angle is in radians
    //center is a [2] array
    // we'll leverage the translate function to make the center the origin
    //lets throw error if points or centers are not arrays
    if (!Array.isArray(points) || !Array.isArray(center)){
        throw new Error(`Invalid type of points or center. Expected Array, got ${typeof points} and ${typeof center}`);
    }

    if (angle_unit == "deg"){
        angle = math.unit(angle, "deg").toNumber("rad");
    }
    // lets throw an error on dimensions mismatch of points and center
    if (points[0].length != 2 || center.length != 2){
        throw new Error("Invalid dimensions of points or center");
    }
    
    let p_T = translate(points, math.multiply(center, -1));
    let R_theta = math.matrix([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]);
    let p_R = math.multiply(p_T, R_theta);
    let p_R_T = translate(p_R.toArray(), center);
    return p_R_T;
}

function scale(points, scaleX, scaleY){
    //points is a [npoints, 2] array, 
    // We must center the points first scale and then translate back to avoid translation when scaling
    // lets throw an error if points is not an array
    if (!Array.isArray(points)){
        throw new Error(`Invalid type of points. Expected Array, got ${typeof points}`);

    }

    let center = math.mean(points, 0);
    let p_T = translate(points, math.multiply(center, -1));
    let S = math.matrix([[scaleX, 0], [0, scaleY]]);
    let p_S = math.multiply(p_T, S);
    let p_S_T = translate(p_S.toArray(), center);
    return p_S_T;
}

function rotateCOM(points, angle, angle_unit="deg"){
    //points is a [npoints, 2] array
    //angle is in radians
    // lets throw an error on dimensions mismatch of points and center
    if (points[0].length != 2){
        throw new Error("Invalid dimensions of points");
    }
    let center = math.mean(points, 0);
    return rotate(points, angle, center, angle_unit);
}



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


function calculateCOMInertiaMoment(points, masses){
    //points is a [npoints, 2] array
    //masses is a [npoints] array
    // lets throw an error on dimensions mismatch of points and masses
    // it is 2d so things can only rotate around the z-axis so we only need the moment of inertia around the z-axis
    // I = sum_i m_i (x_i^2 + y_i^2) = sum_i m_i r_i^2  with r_i = point_i - COM
    if (points.length != masses.length){
        throw new Error("Invalid dimensions of points or masses");
    }
    let COM = calculateCOM(points, masses);
    let pointsCOM = translate(points, math.multiply(COM, -1));
    let r2 = math.sum(pointsCOM.map((p,i) => math.dot(p,p)*masses[i]));
    return r2;
}

function calculateInertiaMoment(points,masses, center){
    if (points.length != masses.length){
        throw new Error("Invalid dimensions of points or masses");
    }
    let pointsCenter = translate(points, math.multiply(center, -1));
    let r2 = math.sum(pointsCenter.map((p,i) => math.dot(p,p)*masses[i]));
    return r2;
}



function _calculateKineticEnergy(vs,masses){
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



//lets write a more general function that will be fed a [ntimesteps,nparticles,2] array of velocities 
//and a [nparticles] array of masses and will return [ntimesteps] array of kinetic energies
// lets make it so if  a [nparticles,2] array of velocities is given it 
//returns a single kinetic energy

function calculateKineticEnergy(vsArray,masses){
    /*
    vsArray is a [ntimesteps,nparticles,2] or [nparticles,2] array of velocities
    masses is a [nparticles] array of masses
    Returns a [ntimesteps] array of kinetic energies or a single kinetic energy if a single timestep is given
    */
    let ndims = getArrayDimensions(vsArray);
    if (ndims == 2){
        //single timestep
        let Ekin = _calculateKineticEnergy(vsArray, masses);
        return Ekin;
    }
    else if (ndims == 3){
        //multiple timesteps
        let EkinArray = vsArray.map(vs => _calculateKineticEnergy(vs, masses));
        return EkinArray;
    }
    else{
        throw new Error(`Invalid dimensions of vsArray. Expected 2 or 3, got ${ndims}`);
    }
}


function getNeighborsDict(STATE){
    let constraintsDistance = STATE.constraints_distance;
    let constraintsPin = STATE.constraints_pin;
    let neighborsDict = {};
    let pinDict = {};
    for(let i = 0; i < constraintsDistance.length; i++){
        let [iParticle, jParticle, distance] = constraintsDistance[i];
        if(neighborsDict[iParticle] == undefined){
            neighborsDict[iParticle] = [[jParticle, distance]];
        }else{
            neighborsDict[iParticle].push([jParticle, distance]);
        }

        if(neighborsDict[jParticle] == undefined){
            neighborsDict[jParticle] = [[iParticle, distance]];
        }else{
            neighborsDict[jParticle].push([iParticle, distance]);
        }
    }

    for(let i = 0; i < constraintsPin.length; i++){
        let [iParticle, pinPoint, distance] = constraintsPin[i];
        pinDict[iParticle] = [pinPoint, distance];
    }

    return [neighborsDict, pinDict];
}



//lets write a couple of functions that will be handy to transform
// model/world coordinates to canvas (centered) coordinates
// basically a shift and a scale, with the y-axis flipped




function model2Canvas(pointsModel, originModel, canvasX, canvasY, modelX, modelY){
    //pointsModel is a [npoints, 2] array
    //originModel is a [2] array, usually [0, 0] but we might want to reference shifted axis in the model space
    //something useful for instance when we want to transform different bodies COM-centered coordinates to canvas coordinates

    //p_{c} = T_{m->c}(p_{m} - o_{m})+o_{c}  where m is model and c is canvas
    //p_{m} = T_{c->m}(p_{c} - o_{c})+o_{m} and T_{c->m} = T_{m->c}^{-1}
    //T_{m->c} = [[sx, 0], [0, -sy]] where sx = canvasX/2 and sy = canvasY/2
    //T_{c->m} = [[1/sx, 0], [0, -1/sy]]
    //lets throw an error if pointsModel or originModel are not arrays
    if (!Array.isArray(pointsModel) || !Array.isArray(originModel)){
        throw new Error(`Invalid type of pointsModel or originModel. Expected Array, got ${typeof pointsModel} and ${typeof originModel}`);
    }
    // lets throw an error if ndims pointsModel and originModel are not the same
    if (pointsModel[0].length != originModel.length){
        throw new Error("Invalid dimensions of pointsModel or originModel");
    }

    let [sx, sy] = [ (canvasX/2)/modelX, (canvasY/2)/modelY];
    let Tmc = math.matrix([[sx, 0], [0, -sy]]);
    let oC = [canvasX/2, canvasY/2];
    let oM = originModel;
    let pM_T = math.transpose( translate(pointsModel, math.multiply(oM, -1)) );
    let pT = math.transpose(math.multiply(Tmc, pM_T));
    let pC = translate(pT.toArray(), oC);
    return pC;
}

function canvas2Model(pointsCanvas, originModel, canvasX, canvasY, modelX, modelY){
    //pointsCanvas is a [npoints, 2] array
    //originModel is a [2] array, usually [0, 0] but we might want to reference shifted axis in the model space

    //p_{c} = T_{m->c}(p_{m} - o_{m})+o_{c}  where m is model and c is canvas
    //p_{m} = T_{c->m}(p_{c} - o_{c})+o_{m} and T_{c->m} = T_{m->c}^{-1}
    //T_{m->c} = [[sx, 0], [0, -sy]] where sx = canvasX/2 and sy = canvasY/2
    //T_{c->m} = [[1/sx, 0], [0, -1/sy]]
    //lets throw an error if pointsCanvas or originModel are not arrays
    if (!Array.isArray(pointsCanvas) || !Array.isArray(originModel)){
        throw new Error(`Invalid type of pointsCanvas or originModel. Expected Array, got ${typeof pointsCanvas} and ${typeof originModel}`);
    }
    // lets throw an error if ndims pointsCanvas and originModel are not the same
    if (pointsCanvas[0].length != originModel.length){
        throw new Error("Invalid dimensions of pointsCanvas or originModel");
    }

    let [sx, sy] = [ (canvasX/2)/modelX, (canvasY/2)/modelY];
    let Tcm =  math.matrix([[1/sx, 0], [0, -1/sy]]);
    let oC = [canvasX/2, canvasY/2];
    let oM = originModel;
    let pC_T = math.transpose( translate(pointsCanvas, math.multiply(oC, -1)) );
    let pT = math.transpose(math.multiply(Tcm, pC_T));
    let pM = translate(pT.toArray(), oM);
    
    return pM;
}


function dArraydtFun(array, h){
    //calculates numericalderivative for an array of values, assuming they come as discrete evaluation of a t->R function sampled at fixed intervals h
    //we'll use a centered difference scheme, for the first and last points we'll use a forward and backward difference scheme
    let dims = getArrayDimensions(array);
    //if dims !=2 ([npoints,ndim]) we throw an error
    if (dims.length != 2){
        let givenDims = dims;
        throw new Error(`Invalid dimensions of array. Expected 2, got ${givenDims}`);
    }

    let dArraydt = [];
    let n = array.length;
    for (let i = 0; i < n; i++){
        if (i == 0){
            dArraydt.push((array[1] - array[0])/h);
        }
        else if (i == n-1){
            dArraydt.push((array[n-1] - array[n-2])/h);
        }
        else{
            dArraydt.push((array[i+1] - array[i-1])/(2*h));
        }
    }

    return dArraydt;
}

function dfdxFun(f, x, h=1e-8){
    return (f(x+h) - f(x-h))/(2*h);
}

function d2fdx2Fun(f, x, h=1e-8){
    return (f(x+h) - 2*f(x) + f(x-h))/(h*h);
}



function bisectionSearch(f,a,b, tolf=1e-4, tolx=1e-4, maxIter=50, verbose=true){
    // we are not finding a root, but a minimum of a function, so we have to do bisection using dfdx as fun not f

    let foriginal = f; // for clarity we'll store the original function
    let [ao, bo] = [a,b];//we'll store the original values of a and b for later checking
    let [funoriginalao, funoriginalbo] = [foriginal(a), foriginal(b)]; // we'll store the original values of f(a) and f(b) for later checking

    let fun = (x) => dfdxFun(f,x); // we are finding the root of this function (minimum of f)

    let fa = fun(a);
    let fb = fun(b);

    // lets add some modifications to ensure we dont find a maximum ( grad = 0 as well)
    // we'll store the minimum value of foriginal found so far, together with x value for that minimum
    let [minf, minx] = [Math.min(funoriginalao, funoriginalbo), funoriginalao < funoriginalbo ? a : b];


    if (verbose) console.log(`a: ${a}, b: ${b}, fa: ${fa}, fb: ${fb}`);

    if(fa*fb >= 0){
        if (verbose) console.log("Root not bracketed returning min(f(a),f(b))");
        return minx;

    }

    let c = (a+b)/2;
    let fc = fun(c);

    let iter = 0;

    while(Math.abs(fc) > tolf && Math.abs(b-a) > tolx && iter < maxIter){
        if (verbose) console.log(`ITER ${iter},   a: ${a}, b: ${b}, c: ${c}, fa: ${fa}, fb: ${fb}, fc: ${fc}`);
        if (fa*fc < 0){
            b = c;
            fb = fc;
        }else{
            a = c;
            fa = fc;
        }

        c = (a+b)/2;
        fc = fun(c);
        iter++;

        // lets check if we have found a new minimum and store it
        if (foriginal(c) < minf){
            minf = foriginal(c);
            minx = c;
        }
    }

    if (iter == maxIter){
        if (verbose) console.log("CONVERGENCE FAILED Max iterations reached");
        //lets throw an error
        throw new Error("CONVERGENCE FAILED Max iterations reached");
    }
    else{
        //lets write a string to the console to show the results
        // we display (xao,foriginal(xao)) and (xbo,foriginal(xbo))  and (minx,foriginal(minx)) together with funder(c) and the number of iterations
        if (verbose) console.log(`a: ${ao}, b: ${bo}, fa: ${funoriginalao},
                                             fb: ${funoriginalbo}, minx: ${minx},
                                              minf: ${minf}, derf(c): ${fc},c: ${c}, iterations: ${iter}`);
        // lets check if the grad root we found correspond to a minimum or a maximum
        // if the second derivative is positive at the root, then we have a minimum else we have a maximum
        if (d2fdx2Fun(f,minx) > 0){
            if (verbose) console.log("Minimum found");
        }else{
            if (verbose) console.log("MAXIMUM FOUND, return minimum value found so far");
        }
        
        return minx;
    }


}


function goldenSectionSearch(f,a,b,tolxrel= 1e-2, maxIter=20, verbose=false){
    
    let phi = (1+Math.sqrt(5))/2;
    let invphi = 1/phi;


    // while (b-a) > tolx
    let c = b - (b-a)*invphi;
    let d = a + (b-a)*invphi;
    let [fa,fc,fd,fb] = [f(a),f(c),f(d),f(b)];
    // lets store the f,x for the minimum of the 4
    let width = b-a;
    let minf = Math.min(fa,fc,fd,fb);
    let minx = minf == fa ? a : minf == fc ? c : minf == fd ? d : b;

    let iter = 0;
    while(Math.abs(b-a)/width > tolxrel && iter < maxIter){
        if (verbose) console.log(`ITER ${iter},   a: ${a}, b: ${b}, c: ${c}, d: ${d}, fc: ${fc}, fd: ${fd}`);
        console.log("Math.abs(b-a)/width",Math.abs(b-a)/width, "tolxrel",tolxrel);
        if (fc < fd){
            b = d;
            fb=fd;

        }
        else{
            a = c;
            fa = fc;
            
        }

        c = b - (b-a)*invphi;
        d = a + (b-a)*invphi;
        fc = f(c);
        fd = f(d);

        let minfn = Math.min(fa,fc,fd,fb);
        let minxn = minf == fa ? a : minf == fc ? c : minf == fd ? d : b;

        if (minfn < minf){
            minf = minfn;
            minx = minxn;
        }

        iter++;
    }

    if (iter == maxIter){
    
        //throw new Error("CONVERGENCE FAILED Max iterations reached");
        console.log("CONVERGENCE FAILED Max iterations reached");
        return minx;
    }

    else{
        if (verbose) console.log(`a: ${a}, b: ${b}, c: ${c}, d: ${d}, fa: ${fa}, fc: ${fc}, fd: ${fd}, fb: ${fb}, minx: ${minx}, minf: ${minf}`);
        return minx;
    }

}



function getNeighborsDelauney(points){
    /*
    Needs the library Delaunator to be loaded on the page
   https://unpkg.com/delaunator@5.0.0/delaunator.min.js
    Example usage
    let neighbors = getNeighbors(points);
    example output
    {
        0: [1, 2, 3],
        1: [0, 2, 3],
        2: [0, 1, 3],
        3: [0, 1, 2]
    }
    */

    const delaunay = Delaunator.from(points);
    const triangles = delaunay.triangles;

    let neighbors = {};
    for (let i = 0; i < triangles.length; i += 3) {
        const p1 = triangles[i];
        const p2 = triangles[i + 1];
        const p3 = triangles[i + 2];
        if (neighbors[p1] === undefined){
            neighbors[p1] = [];
        }
        if (neighbors[p2] === undefined){
            neighbors[p2] = [];
        }
        if (neighbors[p3] === undefined){
            neighbors[p3] = [];
        }
        neighbors[p1].push(p2);
        neighbors[p1].push(p3);
        neighbors[p2].push(p1);
        neighbors[p2].push(p3);
        neighbors[p3].push(p1);
        neighbors[p3].push(p2);
    }

    //lets remove duplicates
    for (let key in neighbors){
        neighbors[key] = [...new Set(neighbors[key])];
    }

    //lets transform keys and values to numbers
    let neighborsNum = {};
    for (let key in neighbors){
        neighborsNum[parseInt(key)] = neighbors[key].map(x => parseInt(x));
    }

    return neighborsNum
}



function lerp(tEval,tArray,xArray){
    //if tEval is outside the range
    // of tArray we throw an error
    let ndim = xArray[0].length;
    let xEval = [];
    //lets print for debuggint teval 0 and final and tarray 0 and final values
    console.log("tEval[0]",tEval[0],"tEval[tEval.length-1]",tEval[tEval.length-1],"tArray[0]",tArray[0],"tArray[tArray.length-1]",tArray[tArray.length-1]);
    for (let i = 0; i < tEval.length; i++){
        let t = tEval[i];
       // if (t < tArray[0] || t > tArray[tArray.length-1]){
         //   throw new Error("tEval out of range");
       // }
        let idx = null;
        for (let j = 0; j < tArray.length; j++){
            if (tArray[j] > t){
                idx = j;
                break;
            }
        }
        //sloppy temporal fix
        if (idx == null){
            idx = tArray.length-2;
        }

        let t1 = tArray[idx-1];
        let t2 = tArray[idx];
        let x1 = xArray[idx-1];
        let x2 = xArray[idx];
        //debuging log
        let x = x1.map((x1i, i) => x1i + (t-t1)*(x2[i]-x1i)/(t2-t1));
        
        xEval.push(x);
    }

    return xEval;
    
}

function EulerStep(fun, x, h,t){
    //fun is the dxdt function fun(x,t)
    //x is a R array
    //h is the step size 

    let dxdt = fun(x,t);
    let xnew = math.add(x, math.multiply(dxdt, h));
    return xnew;
}

function RK4Step(fun, x, h,t){
    //x is a R array
    //h is the step size 

    let k1 = fun(x,t);
    let k2 = fun(math.add(x, math.multiply(k1, h/2)),t+h/2);
    let k3 = fun(math.add(x, math.multiply(k2, h/2)),t+h/2);
    let k4 = fun(math.add(x, math.multiply(k3, h)),t+h);
    let xnew = math.add(x, math.multiply(math.add(k1, math.multiply(k2, 2), math.multiply(k3, 2), k4), h/6));
    return xnew;
}

function RK2Step(fun, x, h, t){
    //x is a R array
    //h is the step size

    let k1 = fun(x,t);
    let k2 = fun(math.add(x, math.multiply(k1, h)),t+h);
    let xnew = math.add(x, math.multiply(math.add(k1, k2), h/2));
    return xnew;
}

function integrateOde(fun, x0, tf, h, method="RK4", tArrayEval = null){
    //fun is the dxdt function, a R x R function where R is the dimension of x
    //x0 is a R array
    //t0 is the initial time
    //tf is the final time
    //h is the step size
    //if tArrayEval is given we use lerping to interpolate the solution at those times
    //method is the integration method, Euler or RK4
    let t0 = 0;
    let x = x0;
    let t = t0;
    let nsteps = Math.floor((tf-t0)/h);
    let xtraj = [x0];
    let ttraj = [t0];
    for (let i = 0; i < nsteps; i++){
        if (method == "Euler"){
            x = EulerStep(fun, x, h, t);
        }
        else if (method == "RK4"){
            x = RK4Step(fun, x, h, t);
        }
        else if (method == "RK2"){
            x = RK2Step(fun, x, h, t);
        }
        else{
            throw new Error("Invalid method, use Euler or RK4");
        }
        t = t + h;
        xtraj.push(x);
        ttraj.push(t);
    }
    if (tArrayEval == null){
        return [xtraj, ttraj];
    }
    else{
        let xEval = lerp(tArrayEval, ttraj, xtraj);
        return xEval;
    }
}

function integrateOdeOnTarray(fun, x0, tArray, method = "RK4"){
    //a simple util wrapper on integrate ODE
    //assumes h is the same for all steps
    let tf = tArray[tArray.length-1];
    let h = tArray[1] - tArray[0];
    let xEval = integrateOde(fun, x0, tf, h, method, tArray);
    return xEval;
}





export {getAngle,translate, rotate, rotateCOM,scale,  model2Canvas, canvas2Model,
     bisectionSearch,goldenSectionSearch, dfdxFun, d2fdx2Fun, dArraydtFun,
      getNeighborsDelauney, calculateCOM, getArrayDimensions, getArrayShape, isShapeEqual, lerp,
      calculateCOMInertiaMoment, calculateInertiaMoment, integrateOde, integrateOdeOnTarray,EulerStep, RK4Step,
    calculateKineticEnergy,getNeighborsDict};