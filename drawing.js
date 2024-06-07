

function drawParticles(svg, xarray, rarray,id_particle, color = "red"){
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
  let particles = svg.selectAll(".particle_"+id_particle).data(data);

  particles.enter().append("circle")
  .attr("class", "particle_"+id_particle)
  .attr("cx", d => d.x)
  .attr("cy", d => d.y)
  .attr("r", d => d.r)
  .attr("id", d => "particle-"+d.idx)
  .style("fill", color);

  particles.attr("cx", d => d.x)
  .attr("cy", d => d.y)
  .attr("r", d => d.r)
  .style("fill", color);


  particles.exit().remove();

  let text = svg.selectAll(".particle_text_"+id_particle).data(data);

  text.enter().append("text")
  .attr("class", "particle_text_"+id_particle)
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
    // Define the arrowhead marker
    // lets group line and head into a svg group
    // lets remove the arrow if it already exists
    svg.append("defs").append("marker")
    .attr("id", arrowId)
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 5)
    .attr("refY", 0)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0,-5L10,0L0,5")
    .attr("class", class_name)
    .attr("fill", color);
    svg.select(`#${arrowId}_line`).remove();



    // Draw a line with the arrowhead marker
    svg.append("line").attr("id", `${arrowId}_line`)
        .attr("x1", x1)
        .attr("y1", y1)
        .attr("x2", x2)
        .attr("y2", y2)
        .attr("stroke-width", 2)
        .attr("stroke", color)
        .attr("marker-end", `url(#${arrowId})`);
        }


function drawState(state,system, svg){
  // draw the particles

  let config = system.getConfig();
  drawParticles(svg, state["xs"], config.masses, "particle", "red");
  // draw the distance constraints
  let constraints_info = system.getConstraintsInfo()
  //constraints_info is an object that might have 2 entries, constraints_distance and constraints
  // if there is constraints_distance we extract particle_indices,distances from this entry
  // if there is constraints_pin, we extract particle_indices,distances,pin_points from this entry
  // and we call drawDistanceConstraints and drawPinConstraints if they exist
  if (constraints_info.hasOwnProperty("constraints_distance")){
    let cd = constraints_info["constraints_distance"];

    drawDistanceConstraints(svg, state["xs"], cd.particle_indices, cd.distances,
                            "constraint", "black");
  }

  if (constraints_info.hasOwnProperty("constraints_pin")){
    let cp = constraints_info["constraints_pin"];
    drawPinConstraints(svg, state["xs"], cp.particle_indices, cp.pin_points, cp.distances, "constraint", "black");
  }

  
  for (let i = 0; i<state["external_forces"].length;i++){
    let [x1,y1] = state["xs"][i];
    let x2 = x1 + state["external_forces"][i][0];
    let y2 = y1 + state["external_forces"][i][1];
    let arrowid = "external_forces_" + i;
    let class_name = "arrow_external_forces";
    drawArrow(svg,x1,y1,x2,y2,arrowid,class_name, "blue");
  }

  for (let i = 0; i<state["constraint_forces"].length;i++){
    let [x1,y1] = state["xs"][i];
    let x2 = x1 + state["constraint_forces"][i][0];
    let y2 = y1 + state["constraint_forces"][i][1];
    let arrowid = "constraint_forces_" + i;
    let class_name = "arrow_constraint_forces";
    drawArrow(svg,x1,y1,x2,y2,arrowid,class_name, "yellow");
  }
}



// lets make a Plot class to handle 1d line plots
  // the idea is that the axis are lazy created, when the first plot is called
  class Plot {
    constructor(div_id, axis_id, x_label = "", y_label = "") {
      this.div_id = div_id;
      this.axis_id = axis_id;
      this.svg = null;
      this.xScale = null;
      this.yScale = null;
      this.x_label = x_label;
      this.y_label = y_label;
      this.legendPosition = 0;// we'll increment it for each plot
   
    }
    createAxis(div_id, axis_id, x_domain, y_domain, x_label = "", y_label=""){
      this.svg = d3.select(`#${div_id}`)
          .append('svg')
          .attr('width', 500)
          .attr('height', 500)
          .attr('id', axis_id);
  
        this.xScale = d3.scaleLinear()
          .domain(x_domain)
          .range([50, 450]);
  
        this.yScale = d3.scaleLinear()
          .domain(y_domain)
          .range([450, 50]);
  
        this.svg.append('g')
          .attr('transform', 'translate(0, 450)')
          .call(d3.axisBottom(this.xScale));
  
        this.svg.append('g')
          .attr('transform', 'translate(50, 0)')
          .call(d3.axisLeft(this.yScale));
  
        this.svg.append('text')
          .attr('x', 250)
          .attr('y', 490)
          .text(x_label)
          .style('text-anchor', 'middle');
  
        this.svg.append('text')
          .attr('x', -250)
          .attr('y', 20)
          .text(y_label)
          .attr('transform', 'rotate(-90)')
          .style('text-anchor', 'middle');
  
  
      }
    plotLine(xarray, yarray, color = "black", label = "") {
      // lets first check if svg is null, if it is not, we plot, if it is, we create axis
      if (this.svg === null) {
        //lets get the x and y domain with extent
        let xdomain = d3.extent(xarray);
        let ydomain = d3.extent(yarray);
        this.createAxis(this.div_id, this.axis_id, 
                       xdomain, ydomain, this.x_label, this.y_label);
      
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
      .attr("x", 460)
      .attr("y", 45 + this.legendPosition)
      .attr("width", 10)
      .attr("height", 10)
      .attr("fill", color);
  
  
      this.svg.append("text")
      .attr("x", 500)
      .attr("y", 50 + this.legendPosition)
      .text(label)
      .style("text-anchor", "middle")
      .style("font-size", "12px")
      .attr("fill", color);
  
      this.legendPosition += 20;
  
  
  
  
  
  
    }
  }
  

  class Drawer{
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
export { drawArrow, drawParticles, drawDistanceConstraints, drawPinConstraints, drawState, drawCross, Plot, Drawer};