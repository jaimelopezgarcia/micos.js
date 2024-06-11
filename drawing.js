


const DEBUG = true;


function drawParticlesOld(svg, xarray, rarray,id_particle,
                           color = "red"){
  // xarray is a [Nparticles, 2] array, rarray is a [Nparticles] array, lets throw an error if they have different lengths
  if (xarray.length != rarray.length){
    throw new Error("xarray and rarray must have the same length");
  }
  // Lets create a data object {x: xarray_x,y: xarray_y, r: rarray}
  // This will be used to create the circles with the enter append exit pattern in d3
  let data = [];
  for (let i = 0; i < xarray.length; i++){
    data.push({x: xarray[i][0], y: xarray[i][1], r: rarray[i],idx: i});
  }
  // lets add the class particle +"_id_particle" and we do the selectall over this class
  // to avoid selecting all the circles in the svg
  let particles = svg.selectAll(".particle"+"_"+id_particle).data(data);

  particles.enter().append("circle")
  .attr("class", "particle_circle")
  .attr("cx", d => d.x)
  .attr("cy", d => d.y)
  .attr("r", d => d.r)
  .attr("id", d => "particle"+"_"+d.idx)
  .style("fill", color);

  particles.attr("cx", d => d.x)
  .attr("cy", d => d.y)
  .attr("r", d => d.r)
  .style("fill", color);


  particles.exit().remove();

  let text = svg.selectAll(".particle_text_"+id_particle).data(data);

  text.enter().append("text")
  .attr("class", "particle_text")
  .attr("x", d => d.x)
  .attr("y", d => d.y)
  .attr("text-anchor", "middle")
  .attr("dominant-baseline", "middle")
  .text(d => d.idx);

  text.attr("x", d => d.x)
  .attr("y", d => d.y);

  text.exit().remove();



  return particles;

}

// Lets write again drawParticles but with several improvements
// We'll get  rid of svg.selectAll and instead we create a group for the circles and text

function drawParticles(svgd3, xs, radii, groupParticlesId,color = "red") {
    // Create or select the group for all particles
    let particlesGroup = svgd3.selectAll(`#${groupParticlesId}`)
        .data([null]) // dummy data
        .join("g")
        .attr("id", groupParticlesId);
      let data = [];
        for (let i = 0; i < xs.length; i++){
          data.push({x: xs[i][0], y: xs[i][1], r: radii[i],idx: i});
        }
    // Bind data to existing particle groups
    let particles = particlesGroup.selectAll(".particle")
        .data(data);

    // Handle new particles (enter)
    let newParticles = particles.enter()
        .append("g")
        .attr("class", "particle")
        .attr("id", d => `${groupParticlesId}_particle_${d.idx}`);

    newParticles.append("circle");
    newParticles.append("text");

    particles = newParticles.merge(particles);



    // Update existing particles (update)
    particles.select("circle")
        .attr("cx", (d, i) => d.x)
        .attr("cy", (d, i) => d.y)
        .attr("r", (d, i) => d.r)
        .attr("fill", color);

    particles.select("text")
        .attr("x", (d, i) => d.x)
        .attr("y", (d, i) => d.y)
        .text((d, i) => `${i}`);

    // Remove old particles (exit)
    particles.exit().remove();
  }



// This function will draw the constraints between the particles
function drawDistanceConstraints(svg, xarray, constraints_idxs,constraints_distances,
                                 id_constraint, color = "black"){
// xarray is a [Nparticles,2] array, constraints_idxs is a [Nconstraints,2] array with the indices
// of the particles that are connected by the constraint, and
// constraints_distances is a [Nconstraints] array of the distances that must be enforced
// if constraint is satisfied the distance between the particles will
// be equal to the distance in constraints_distances

// Lets create a data object 
//{x1: x1, y1: y1, x2: x2, y2: y2, distance: distance,ratios: ratios}
// lets adopt the convention that x2-x1 where x1 is the first particle in the constraint
let data = [];
for (let i = 0; i < constraints_idxs.length; i++){
  let idxs = constraints_idxs[i];

  let x1 = xarray[idxs[0]][0];
  let y1 = xarray[idxs[0]][1];
  let x2 = xarray[idxs[1]][0];
  let y2 = xarray[idxs[1]][1];
  let distance = constraints_distances[i];
  let realdistance = Math.sqrt((x2-x1)**2 + (y2-y1)**2);
  // lets draw the distance constraint with the imposed distance instead of the real distance
  // lets calculate the \vec{r} = (x2-x1,y2-y1) and the
  let runit_x = (x2-x1)/realdistance;
  let runit_y = (y2-y1)/realdistance;
  let x2_imposed = x1 + runit_x*distance;
  let y2_imposed = y1 + runit_y*distance;
  data.push({x1: x1, y1: y1, x2: x2_imposed, y2: y2_imposed, distance: distance});

}
// lets add the class constraint +"_id_constraint" and we do the selectall over this class
// to avoid selecting all the lines in the svg
let constraints = svg.selectAll(".distance_constraint_"+id_constraint).data(data);

constraints.enter().append("line")
.attr("class", "distance_constraint_"+id_constraint)
.attr("x1", d => d.x1)
.attr("y1", d => d.y1)
.attr("x2", d => d.x2)
.attr("y2", d => d.y2)
.attr("stroke", color);

constraints.attr("x1", d => d.x1)
.attr("y1", d => d.y1)
.attr("x2", d => d.x2)
.attr("y2", d => d.y2)
.attr("stroke", color);

constraints.exit().remove();



return constraints;
                                 }

                                
function drawCross(svg, x, y, size, id, color = "black") {
  // lets remove the cross if it already exists
    svg.select(`#${id}`).remove();
    var cross = svg.append("g").attr("id",id)
    cross.append("line")
    .attr("x1", x - size / 2)
    .attr("y1", y - size / 2)
    .attr("x2", x + size / 2)
    .attr("y2", y + size / 2)
    .attr("stroke", color)


    cross.append("line")
    .attr("x1", x - size / 2)
    .attr("y1", y + size / 2)
    .attr("x2", x + size / 2)
    .attr("y2", y - size / 2)
    .attr("stroke", color)

}
function drawPinConstraints(svg, xarray, constraints_idxs,
                            pinpoints ,distances, id_constraint, color = "black"){
  //very similar to drawDistanceConstraints but now there is only 1 particle involved
  // xarray [Nparticles,2] array, constraints_idxs [Nconstraints] pinpoints is a [Nconstraints,2] 
  //the other coordinates correspond to the pin point
  // Lets throw an error if different lengths of pinpoints and constraints_idxs
  if (constraints_idxs.length != pinpoints.length){
    throw new Error("constraints_idxs and pinpoints must have the same length");
  }
  // Lets create a data object {p1:p1,p2:p2,x1:x1,x2:x2, distance: distance,ratios: ratios}
  // Lets adopt the convention that x-p where x is the particle and p is the pin point
  let data = [];
  for (let i = 0; i < constraints_idxs.length; i++){
    let idx = constraints_idxs[i];
    let x1 = xarray[idx][0];
    let y1 = xarray[idx][1];
    let x2 = pinpoints[i][0];
    let y2 = pinpoints[i][1];
    let distance = distances[i];
    let realdistance = Math.sqrt((x2-x1)**2 + (y2-y1)**2);
    // lets draw the distance constraint with the imposed distance instead of the real distance
    // lets calculate the \vec{r} = (x2-x1,y2-y1) and the
    let runit_x = (x2-x1)/realdistance;
    let runit_y = (y2-y1)/realdistance;
    let x2_imposed = x1 + runit_x*distance;
    let y2_imposed = y1 + runit_y*distance;
    data.push({x1: x1, y1: y1, x2i: x2_imposed, y2i: y2_imposed,
              x2:x2,y2:y2, distance: distance});
  }
  // lets add the class pin_constraint +"_id_constraint" and we do the selectall over this class
  // to avoid selecting all the lines in the svg
  let constraints = svg.selectAll(".pin_constraint_"+id_constraint).data(data);

  constraints.enter().append("line")
  .attr("class", "pin_constraint_"+id_constraint)
  .attr("x1", d => d.x1)
  .attr("y1", d => d.y1)
  .attr("x2", d => d.x2i)
  .attr("y2", d => d.y2i)
  .attr("stroke", color);

  constraints.attr("x1", d => d.x1)
  .attr("y1", d => d.y1)
  .attr("x2", d => d.x2i)
  .attr("y2", d => d.y2i)
  .attr("stroke", color);

  constraints.exit().remove();

  // additionally lets draw pin point with a cross, its coordinates are x2,y2
  // using the function drawCross
  for (let i = 0; i<data.length;i++){
    let x2 = data[i].x2;
    let y2 = data[i].y2;
    let id_cross = "pin_cross_"+i;
    drawCross(svg, x2, y2, 10, id_cross, color = "black");
  }




  return constraints;

}





function drawArrow(svg, x1, y1, x2, y2, arrowId,class_name = "arrow", color="black") {

    // lets remove the arrow if it already exists
    svg.select(`#${arrowId}`).remove();
    // lets group line and head into a svg group
    let arrow = svg.append("g").attr("id",arrowId);
    // Draw a line with the arrowhead marker
    arrow.append("line")
        .attr("x1", x1)
        .attr("y1", y1)
        .attr("x2", x2)
        .attr("y2", y2)
        .attr("stroke-width", 2)
        .attr("stroke", color)
        .attr("class", class_name)
        .attr("id", `${arrowId}_line`);
    // Define the arrowhead marker
    arrow.append("defs").append("marker")
    .attr("id", `${arrowId}_head`)
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 5)
    .attr("refY", 0)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0,-5L10,0L0,5")
    .attr("class", class_name)
    .attr("fill", color)
    .attr("id", `${arrowId}_head`);
    // lets add the marker to the line
    arrow.select("line").attr("marker-end", `url(#${arrowId}_head)`);

    }






  class Plot{
    constructor(div_id, axis_id,
       x_label = "", y_label = "", width = 500, height = 500) {
      this.div_id = div_id;
      this.axis_id = axis_id;
      this.svg = null; // lazy creation when first plot is called
      this.xScale = null;
      this.yScale = null;
      this.x_label = x_label;
      this.y_label = y_label;
      this.legendPosition = 0;// we'll increment it for each plot
      this.width = width;
      this.height = height;
      //lets create a list of colors to cycle through if color is not provided
      this.colors = ["black", "blue", "red", "green", "purple", "orange", "brown", "gray", "pink", "cyan", "magenta"];
      this.counter = 0;
    }
    createAxis(div_id, axis_id, x_domain, y_domain, x_label = "", y_label=""){

      
      this.svg = d3.select(`#${div_id}`)
          .append('svg')
          .attr('width', this.width)
          .attr('height', this.height)
          .attr('id', axis_id);

    

      

      this.xScale = d3.scaleLinear()
        .domain(x_domain)
        .range([50, this.width - 50]);

      this.yScale = d3.scaleLinear()
        .domain(y_domain)
        .range([this.height - 50, 50]);

      this.svg.append('g')
        .attr('transform', `translate(0, ${this.height - 50})`)
        .call(d3.axisBottom(this.xScale));

      this.svg.append('g')
        .attr('transform', `translate(50, 0)`)
        .call(d3.axisLeft(this.yScale));

      this.svg.append('text')
        .attr('x', this.width / 2)
        .attr('y', this.height - 10)
        .text(x_label)
        .style('text-anchor', 'middle');

      this.svg.append('text')
        .attr('x', -this.height / 2)
        .attr('y', 20)
        .text(y_label)
        .attr('transform', 'rotate(-90)')
        .style('text-anchor', 'middle');
    
  

    }
    plotLine(xarray, yarray,label = "", color = null) {
      // lets first check if svg is null, if it is not, we plot, if it is, we create axis
      if (this.svg === null) {
        //lets get the x and y domain with extent
        if (DEBUG){
          console.log("Plotting for the first time, creating axis for div_id", this.div_id)
        }
        let xdomain = d3.extent(xarray);
        let ydomain = d3.extent(yarray);
        this.createAxis(this.div_id, this.axis_id, 
                       xdomain, ydomain, this.x_label, this.y_label);
      
      }
      if (color === null){
        color = this.colors[this.counter];
        this.counter = (this.counter + 1) % this.colors.length;
        if (DEBUG){
          console.log("Color not provided, using color", color)
        }
      }

  
      let line = d3.line()
      .x(d => this.xScale(d.x))
      .y(d => this.yScale(d.y));
  
      let data = xarray.map((x,i) => {return {x:x, y:yarray[i]} } );
  
      // lets build the id with axis_id and label
      let id_data  = `${this.axis_id}-${label}`
      this.svg.append("path")
        .datum(data)
        .attr("fill", "none")
        .attr("stroke", color)
        .attr("stroke-width", 1.5)
        .attr("d", line)
        .attr("id", id_data);
  
      // lets add the label to the legend
  
      this.svg.append("rect")
      .attr("x", this.width - 40)
      .attr("y", 45 + this.legendPosition)
      .attr("width", 10)
      .attr("height", 10)
      .attr("fill", color);
  
  
      this.svg.append("text")
      .attr("x", this.width - 20)
      .attr("y", 50 + this.legendPosition)
      .text(label)
      .style("text-anchor", "middle")
      .style("font-size", "12px")
      .attr("fill", color);
  
      this.legendPosition += 20;
  
  
  
  
    }
    clearPlot(){
      this.svg.selectAll("*").remove();
      this.svg = null;
      this.legendPosition = 0;
    }
    }
  


function drawForces(STATE, svg){

      //lets throw an error if the forces are not in the state and output the missing forces
      let name_forces = ["external_forces", "constraint_forces", "contact_forces", "friction_forces"];
      let missing_forces = name_forces.filter(name => !STATE.hasOwnProperty(name));
      if (missing_forces.length > 0){
          throw new Error("Missing forces in STATE: " + missing_forces);
      }

      // error if state missing xs
      if (!STATE.hasOwnProperty("xs")){
          throw new Error("Missing xs in STATE");
      }

      
      let fext = STATE.external_forces;
      let fconstraint = STATE.constraint_forces;
      let fcontact = STATE.contact_forces;
      let ffriction = STATE.friction_forces;
      let xs = STATE.xs;
      for (let i = 0; i < xs.length; i++){
          let [x,y] = xs[i];
          let [fx, fy] = fext[i];
          drawArrow(svg, x, y, x+fx, y+fy, "arrow_fext_" + i, "arrow fext", "blue");
          [fx, fy] = fconstraint[i];
          drawArrow(svg, x, y, x+fx, y+fy, "arrow_fconstraint_" + i, "arrow fconstraint", "green");
          [fx, fy] = fcontact[i];
          drawArrow(svg, x, y, x+fx, y+fy, "arrow_fcontact_" + i, "arrow fcontact", "red");
          [fx, fy] = ffriction[i];
          drawArrow(svg, x, y, x+fx, y+fy, "arrow_ffriction_" + i, "arrow ffriction", "purple");
      }
  }


  function drawConstraints(STATE, svg){
    /*
    *   Function to draw the constraints in the svg
    *   The constraints are in the STATE object
    *  The constraints are of 2 types, distance and pin
    * The distance constraints are of the form [idx1, idx2, distance]
    * The pin constraints are of the form [idx, [x,y], distance]
    */
    let isConsDistEmpty = STATE.constraints_distance.length === 0;
    let isConsPinEmpty = STATE.constraints_pin.length === 0;
    let xs = STATE.xs;
    if (!isConsDistEmpty){
      let constraints_distance = STATE.constraints_distance;
      let distances1 = constraints_distance.map(constraint => constraint[2]);
      let constraints_idxs1 = constraints_distance.map(constraint => [constraint[0], constraint[1]]);
      drawDistanceConstraints(svg, xs, constraints_idxs1, distances1, "distance_constraints", "black");
    }
    if (!isConsPinEmpty){
      let constraints_pin = STATE.constraints_pin;
      let constraints_idxs2 = constraints_pin.map(constraint => constraint[0]);
      let pinpoints = constraints_pin.map(constraint => constraint[1]);
      let distances2 = constraints_pin.map(constraint => constraint[2]);
      drawPinConstraints(svg, xs, constraints_idxs2, pinpoints, distances2, "pin_constraints", "black");
    }
    
}


function plotStateStoryPositionsVelocities(ParentDivId,STATE_STORY, VISUAL_OPTIONS){

  let parentDiv = document.getElementById(ParentDivId);
  let tarray = STATE_STORY.map(state => state.time);
  
  let nparticles = STATE_STORY[0].xs.length;
  let y_names = ["xplot", "yplot", "vxplot", "vyplot"];
  let div_plots_ids = [];

  for (let y_name of y_names){
      let plot_id = ParentDivId + "-" + y_name;
      let div = parentDiv.querySelector("#" + plot_id);
      if (div==null){
          div = document.createElement("div");
          div.id = plot_id;
          parentDiv.appendChild(div);
      }
      div_plots_ids.push(div.id);

  }
  
  // we need 4 plots, on each plot we'll plot all particles
  let tracesx = [];
  let tracesy = [];
  let tracesvx = [];
  let tracesvy = [];
  for (let i = 0; i < nparticles; i++){
      let trace_x = {
          x: tarray,
          y: STATE_STORY.map(state => state.xs[i][0]),
          mode: "lines",
          name: "particle_" + i + "_x",
      }
      let trace_y = {
          x: tarray,
          y: STATE_STORY.map(state => state.xs[i][1]),
          mode: "lines",
          name: "particle_" + i + "_y",
      }
      let trace_vx = {
          x: tarray,
          y: STATE_STORY.map(state => state.vs[i][0]),
          mode: "lines",
          name: "particle_" + i + "_vx",
      }
      let trace_vy = {
          x: tarray,
          y: STATE_STORY.map(state => state.vs[i][1]),
          mode: "lines",
          name: "particle_" + i + "_vy",
      }
      tracesx.push(trace_x);
      tracesy.push(trace_y);
      tracesvx.push(trace_vx);
      tracesvy.push(trace_vy);
  }

  let layout = {
      title: "State Story",
      xaxis: {
          title: "time",
      },
      yaxis: {
          title: "value",
      },
      width: VISUAL_OPTIONS["width_plots"],
      height: VISUAL_OPTIONS["height_plots"],
      margin: {
          l: 10,
          r: 10,
          b: 20,
          t: 20,
          pad: 4,
      }
  }

  Plotly.newPlot(div_plots_ids[0], tracesx, layout);
  Plotly.newPlot(div_plots_ids[1], tracesy, layout);
  Plotly.newPlot(div_plots_ids[2], tracesvx, layout);
  Plotly.newPlot(div_plots_ids[3], tracesvy, layout);

 
}


  class DrawerOLD{
    constructor(svg, width, height, particle_density){
      this.svg = svg;
      this.width = width;
      this.height = height;
      this.x_model_domain = [-1,1];
      this.y_model_domain = [-1,1];
      this.origin_model = [0,0];
      this.origin_canvas = [width/2,height/2];// canvas point where the origin will be placed
      this.scale_factor_x = width/(this.x_model_domain[1] - this.x_model_domain[0]);
      this.scale_factor_y = -height/(this.y_model_domain[1] - this.y_model_domain[0]);//negative because in svg y grows down
      this.particle_density = particle_density;
    }
    model2Canvas(array, shift = false){
      //lets check type, array is [N,2] model data is expected to be O(1)
      //shift is an optional parameter because coord points will be shifted by the origin but forces will not ( second derivative of the position)
      if (!Array.isArray(array)){
        throw new Error("Array should be an array")
      }
      if (array[0].length === 0){
        throw new Error("Array must be of shape [N,2]")
      }
  
      let array_canvas = array.map(point => {
        let x = point[0];
        let y = point[1];
        let origin_model  = [0,0];
        let origin_canvas = [0,0];
        if (shift){
          origin_model  = [this.origin_model[0],this.origin_model[1]];
          origin_canvas = [this.origin_canvas[0],this.origin_canvas[1]];
        }
        else{
          origin_model  = [0,0];
          origin_canvas = [0,0];
        }
        let x_canvas = (x + origin_model[0]) * this.scale_factor_x + origin_canvas[0];
        let y_canvas = (y + origin_model[1]) * this.scale_factor_y + origin_canvas[1];
        return [x_canvas,y_canvas];
      } );
  
      return array_canvas;
    }
    canvas2Model(array){
      //lets check type, array is [N,2] 
      if (!Array.isArray(array)){
        throw new Error("Array should be an array")
      }
      if (array[0].length === 0){
        throw new Error("Array must be of shape [N,2]")
      }
  
      let array_model = array.map(point => {
        let x_canvas = point[0];
        let y_canvas = point[1];
        let x = (x_canvas + this.origin_canvas[0]) / this.scale_factor_x + this.origin_model[0];
        let y = (y_canvas + this.origin_canvas[1]) / this.scale_factor_y + this.origin_model[1];
        return [x,y];
      });
  
      return array_model;
    
    }
  
    _calculateDistance(point1, point2){
      // constraint distances are calculated in the model domain, we cant turn the distance model scalar
      // into canvas distance, we need the points in the canvas domain
      return Math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2);
  
  
    }
  
    //lets write the different drawing functions that are in drawState2
  
    _drawArrowsArray(svg, xs_origin_array, arrows_array, arrows_id, class_name, color){
  
      for ( let i = 0; i<xs_origin_array.length;i++){
        let [x1,y1] = xs_origin_array[i];
        let x2 = x1 + arrows_array[i][0];
        let y2 = y1 + arrows_array[i][1];
        let arrowid = arrows_id + i;

        drawArrow(svg,x1,y1,x2,y2,arrowid,class_name, color);
      }
    }
  
    _masses2Radiuses(masses, particle_density){
      let radiuses = masses.map(mass=>Math.sqrt(mass/particle_density/Math.PI));
      return radiuses;
    }
    drawState(state, system){
      let xs_model = state["xs"];
      let external_forces_model = state["external_forces"];
      let constraint_forces_model = state["constraint_forces"];
  
      let xs_canvas = this.model2Canvas(xs_model, true);
      let external_forces_canvas = this.model2Canvas(external_forces_model, false);
      let constraint_forces_canvas = this.model2Canvas(constraint_forces_model, false);
  
      let radiuses = this._masses2Radiuses(system.masses, this.particle_density);

      drawParticles(this.svg, xs_canvas, radiuses, "particle", "red");

      this._drawArrowsArray(this.svg, xs_canvas,
                               external_forces_canvas,
                             "external_forces", "arrows_external_forces", "blue");
  
      this._drawArrowsArray(this.svg, xs_canvas,
                               constraint_forces_canvas,
                              "constraint_forces", "arrows_constraint_forces", "yellow");
      
      let constraints_info = system.getConstraintsInfo();
  
      if (constraints_info.hasOwnProperty("constraints_distance")){
          let cd = constraints_info["constraints_distance"];
          let pi = cd.particle_indices;
          let distances_canvas = pi.map((indices,i)=>{
            let [idx1,idx2] = indices;
            let point1 = xs_canvas[idx1];
            let point2 = xs_canvas[idx2];
            let distance = this._calculateDistance(point1,point2);
            return distance;
          });
          drawDistanceConstraints(svg, xs_canvas, pi, distances_canvas,
                               "constraint", "black");
    }
  
    if (constraints_info.hasOwnProperty("constraints_pin")){
      let cp = constraints_info["constraints_pin"];
      let pi = cp.particle_indices;
      let pin_points_canvas = this.model2Canvas(cp.pin_points,true);
      let distances_canvas = pi.map((indices,i)=>{
        let [idx1,idx2] = indices;
        let point1 = xs_canvas[idx1];
        let point2 = pin_points_canvas[i];
        let distance = this._calculateDistance(point1,point2);
        return distance;
      });
      drawPinConstraints(svg,xs_canvas, pi, pin_points_canvas, distances_canvas, "constraint", "black");
    }
                           
  
    }
  }


  class Drawer{
    /*

    Lets make some changes, lets extract width and height from the svg dynamically
    */
    constructor(svg, x_model_domain = [-1,1], y_model_domain = [-1,1], origin_model = [0,0]){
      this.svg = svg;// html dom, not d3
      this.x_model_domain = x_model_domain;
      this.y_model_domain = y_model_domain;
      this.origin_model = origin_model;

    }

    getCanvasDimensions(){
      let width = this.svg.clientWidth;
      let height = this.svg.clientHeight;
      let origin_canvas = [width/2,height/2];// canvas point where the origin will be placed
      let scale_factor_x = width/(this.x_model_domain[1] - this.x_model_domain[0]);
      let scale_factor_y = -height/(this.y_model_domain[1] - this.y_model_domain[0]);//negative because in svg y grows down

      return [width,height,origin_canvas,scale_factor_x,scale_factor_y];
    }
    model2Canvas(array, shift = false){
      //lets check type, array is [N,2] model data is expected to be O(1)
      //shift is an optional parameter because coord points will be shifted by the origin but forces will not ( second derivative of the position)
      if (!Array.isArray(array)){
        throw new Error("Array should be an array")
      }
      if (array[0].length === 0){
        throw new Error("Array must be of shape [N,2]")
      }

      let [width,height,origin_canvas,scale_factor_x,scale_factor_y] = this.getCanvasDimensions();
  
      let array_canvas = array.map(point => {
        let x = point[0];
        let y = point[1];
        let origin_model  = [0,0];
        if (shift){
          origin_model  = [this.origin_model[0],this.origin_model[1]];
          origin_canvas = [origin_canvas[0],origin_canvas[1]];
        }
        else{
          origin_model  = [0,0];
          origin_canvas = [0,0];
        }
        let x_canvas = (x + origin_model[0]) * scale_factor_x + origin_canvas[0];
        let y_canvas = (y + origin_model[1]) * scale_factor_y + origin_canvas[1];
        return [x_canvas,y_canvas];
      } );
  
      return array_canvas;
    }
    canvas2Model(array,shift = false){
      //lets check type, array is [N,2] 
      if (!Array.isArray(array)){
        throw new Error("Array should be an array")
      }
      if (array[0].length === 0){
        throw new Error("Array must be of shape [N,2]")
      }

      let [width,height,origin_canvas,scale_factor_x,scale_factor_y] = this.getCanvasDimensions();
  
      let array_model = array.map(point => {
        let x_canvas = point[0];
        let y_canvas = point[1];
        let origin_model  = [0,0];
        if (shift){
          origin_model  = [this.origin_model[0],this.origin_model[1]];
          origin_canvas = [origin_canvas[0],origin_canvas[1]];
        }
        else{
          origin_model  = [0,0];
          origin_canvas = [0,0];
        }
        let x = (x_canvas - origin_canvas[0]) / scale_factor_x + origin_model[0];
        let y = (y_canvas - origin_canvas[1]) / scale_factor_y + origin_model[1];
        return [x,y];
      });
  
      return array_model;
    
    }
  
    _calculateDistance(point1, point2){
      // constraint distances are calculated in the model domain, we cant turn the distance model scalar
      // into canvas distance, we need the points in the canvas domain
      return Math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2);
  
  
    }
  

  
    _masses2Radiuses(masses, particle_density){
      let radiuses = masses.map(mass=>Math.sqrt(mass/particle_density/Math.PI));
      return radiuses;
    }

    state2Canvas(STATE){
      // Here we transform the different magnitudes in STATE to canvas coordinates, we must generate a twin canvasSTATE object
      //lets start by deep copying the STATE object
      // positions are transformed with  shift = true, forces and velocities are transformed with shift = false

      let canvasSTATE = JSON.parse(JSON.stringify(STATE));
      canvasSTATE.xs = this.model2Canvas(STATE.xs,true);
      canvasSTATE.vs = this.model2Canvas(STATE.vs,false);
      canvasSTATE.external_forces = this.model2Canvas(STATE.external_forces,false);
      canvasSTATE.constraint_forces = this.model2Canvas(STATE.constraint_forces,false);
      canvasSTATE.contact_forces = this.model2Canvas(STATE.contact_forces,false);
      canvasSTATE.friction_forces = this.model2Canvas(STATE.friction_forces,false);

      //we need to transform constraints data, constraints_distance is a list of [idx1,idx2,distance]
      // we need to recalculated the distance in canvas coordinates
      let canvasConstraintsDistance = canvasSTATE.constraints_distance.map(constraint => {
        let [idx1,idx2,distance] = constraint;
        //we need to get in model coordinates [dx,dy] calculate the current distance ( not equal to the distance in the constraint unless the constraint is satisfied)
        //get the point at distance/currendistance * [dx,dy] from point1 and transform it to canvas coordinates
        //once this is done we calculate distance_canvas
        let point1 = STATE.xs[idx1];
        let point2 = STATE.xs[idx2];
        let dx = point2[0] - point1[0];
        let dy = point2[1] - point1[1];
        let current_distance = Math.sqrt(dx**2 + dy**2);
        let ratio = distance/current_distance;
        let canvasdr = this.model2Canvas([[ratio*dx,ratio*dy]],false)[0];
        let distance_canvas = Math.sqrt(canvasdr[0]**2 + canvasdr[1]**2);

        return [idx1,idx2,distance_canvas];
      });

      canvasSTATE.constraints_distance = canvasConstraintsDistance;

      // now the same for the pin constraints points and distances
      let canvasConstraintsPin = canvasSTATE.constraints_pin.map(constraint => {
        let [idx,pinpoint,distance] = constraint;
        // we employ the same strategy, we calculate the distance in model coordinates, then we calculate the point at distance from the particle
        // and we transform it to canvas coordinates
        let point1 = STATE.xs[idx];
        let dx = pinpoint[0] - point1[0];
        let dy = pinpoint[1] - point1[1];
        let current_distance = Math.sqrt(dx**2 + dy**2);
        let ratio = distance/current_distance;
        let canvasdr = this.model2Canvas([[ratio*dx,ratio*dy]],false)[0];
        let distance_canvas = Math.sqrt(canvasdr[0]**2 + canvasdr[1]**2);
        let canvasPinPoint = this.model2Canvas([pinpoint],true)[0];
        

        return [idx,canvasPinPoint,distance_canvas];
      });

      canvasSTATE.constraints_pin = canvasConstraintsPin;

      return canvasSTATE;


    }

    drawState(STATE, VISUAL_OPTIONS, particles_group_id){
      let s = STATE;
      let svgd3 = d3.select(this.svg);
      let particle_density = VISUAL_OPTIONS["particle_density"];
      let radii = s.masses.map(mass => Math.sqrt(mass/particle_density));
      
      let canvasSTATE = this.state2Canvas(s);

      //lets draw particles
      drawParticles(svgd3, canvasSTATE.xs, radii, particles_group_id, "red");
      //lets draw forces
      drawForces(canvasSTATE, svgd3);
      //lets draw constraints
      drawConstraints(canvasSTATE, svgd3);
  

    }
    //lets draw a scale reference in a corner of the svg, a square with a side of 0.1 in model coordinates
    drawScaleReference(){
      //lets write down the square in model coords and then we transform it to canvas coords
      //lets draw it in the top left corner, using relative measure to this.y_model_domain and this.x_model_domain
      let square_model = [[0,0],[0.1,0],[0.1,0.1],[0,0.1],[0,0]];
      //lets shift it in the model coords to the top left corner, top left corner is this.x_model_domain[0], this.y_model_domain[1]
      let square_model_shifted = square_model.map(point => [point[0] + this.x_model_domain[0]+0.2, point[1] + this.y_model_domain[1]-0.2]);
      //lets transform it to canvas coords
      let square_canvas = this.model2Canvas(square_model_shifted,true);
      //lets draw the square
      let svgd3 = d3.select(this.svg);
      svgd3.append("path")
      .datum(square_canvas)
      .attr("fill", "none")
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .attr("d", d3.line()
      .x(d => d[0])
      .y(d => d[1])
      .curve(d3.curveLinearClosed));
      //lets add a text with the scale, lets add some informative text to the right of the square
      svgd3.append("text")
      .attr("x", square_canvas[0][0]+20)
      .attr("y", square_canvas[0][1]+10)
      .attr("fill", "brown")
      .attr("font-size", "12px")
      .attr("text-anchor", "middle")
      .attr("alignment-baseline", "middle")
      .text("0.1 square model units");



    }

  }
export { drawArrow, drawParticles, drawDistanceConstraints,
   drawPinConstraints, drawCross, drawForces,drawConstraints,plotStateStoryPositionsVelocities , Plot, Drawer};