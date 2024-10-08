import {Polygon} from "./solver.js";
import {getSvgRelativeCoords} from "./utils.js";
//lets import external dependencies: d3, mathjs,
import * as d3 from "d3";
import * as math from "mathjs";


const DEBUG = true;


/*Lets write a very succint description of the functions in this file, basically constructor signatures 
and a minimalistic description, lets skip utility private-like function as for instance newSvgElmnt

newSvgElmnt(canvas,tag) -> creates a new svg element with tag and appends it to canvas
drawCircles(svg,xs, radii = 5, color = "red", groupCirclesId = "circlesGroup", display_indices = true) -> draws circles in svg
drawLines(svg, xsOrigs, xsDests, color = "black", groupLinesId = "linesGroup") -> draws lines in svg
drawArrow(svg,porigin,pend, id = "arrow", class_name = "arrow",color = "black", width = 2) -> draws an arrow in svg
drawParticles(svgd3, xs, radii, groupParticlesId,color = "red") -> draws particles in svgd3
drawDistanceConstraints(svg, xarray, constraints_idxs,constraints_distances,id_constraint, color = "black") -> draws distance constraints in svg
drawCross(svg, x, y, size, id, color = "black") -> draws a cross in svg
drawPinConstraints(svg, xarray, constraints_idxs,pinpoints ,distances, id_constraint, color = "black") -> draws pin constraints in svg
drawArrowd3(svgd3, x1, y1, x2, y2, arrowId,class_name = "arrow", color="black") -> draws an arrow in svgd3
drawForces(STATE, svg, scaleArrows = 1.0) -> draws forces in svg
drawConstraints(STATE, svg) -> draws constraints in svg
plotStateStoryPositionsVelocities(ParentDivId,STATE_STORY, VISUAL_OPTIONS) -> plots state story in div
Plot -> class to plot in divs
Drawer -> class to transform model to canvas coordinates

*/

function newSvgElmnt(canvas,tag){
  let element = document.createElementNS('http://www.w3.org/2000/svg', tag);
  canvas.appendChild(element);
  return element;
}
// lets write a function to draw circles, but with plain js not d3

function drawCircles(svg,xs, radii = 5, color = "red",
             groupCirclesId = "circlesGroup", display_indices = true){

    // lets remove the group if it already exists
    let groupCircles = svg.querySelector(`#${groupCirclesId}`);
    if (groupCircles != null){
        groupCircles.remove();
    }
    // lets create the group
    groupCircles = document.createElementNS("http://www.w3.org/2000/svg", "g");
    groupCircles.id = groupCirclesId;
    svg.appendChild(groupCircles);
    for (let i = 0; i < xs.length; i++){
      //lets create  a group for the circle and the text

        let gInner = newSvgElmnt(groupCircles,"g");
        gInner.id = `${groupCirclesId}_${i}`;

        let x = xs[i][0];
        let y = xs[i][1];
        let circle = newSvgElmnt(gInner,"circle");
        circle.setAttribute("cx", x);
        circle.setAttribute("cy", y);
        circle.setAttribute("r", radii);
        circle.setAttribute("fill", color);
        if (display_indices){
            let text = newSvgElmnt(gInner,"text");
            text.setAttribute("x", x);
            text.setAttribute("y", y);
            text.textContent = i;
        }
    }

}
// again, a simple funciton to draw lines with plain js

function drawLines(svg, xsOrigs, xsDests, color = "black", groupLinesId = "linesGroup"){
    // lets remove the group if it already exists
    let groupLines = svg.querySelector(`#${groupLinesId}`);
    if (groupLines != null){
        groupLines.remove();
    }
    // lets create the group
    groupLines = newSvgElmnt(svg,"g");
    groupLines.id = groupLinesId;
    for (let i = 0; i < xsOrigs.length; i++){
        let x1 = xsOrigs[i][0];
        let y1 = xsOrigs[i][1];
        let x2 = xsDests[i][0];
        let y2 = xsDests[i][1];
        let line = newSvgElmnt(groupLines,"line");
        line.setAttribute("x1", x1);
        line.setAttribute("y1", y1);
        line.setAttribute("x2", x2);
        line.setAttribute("y2", y2);
        line.setAttribute("stroke", color);
    }
}


function drawArrow(svg,porigin,pend, id = "arrow", class_name = "arrow",
                     color = "black", width = 2){
            
    //lets create a group for the arrow line and the arrow head
    let group = newSvgElmnt(svg,"g");
    group.setAttribute("id",id)
    group.setAttribute("class",class_name);
    //lets draw the line

    let [xo,yo] = porigin;
    let [xe,ye] = pend;
    let line = newSvgElmnt(group,"line");
    line.setAttribute("x1",xo);
    line.setAttribute("y1",yo);
    line.setAttribute("x2",xe);
    line.setAttribute("y2",ye);
    line.setAttribute("stroke",color);
    line.setAttribute("stroke-width",width);

    //lets draw the arrow head

    let markerHead = newSvgElmnt(svg,"marker");
    markerHead.setAttribute("id","arrowhead");
    markerHead.setAttribute("markerWidth","20");
    markerHead.setAttribute("markerHeight","20");
    markerHead.setAttribute("refX","0");
    markerHead.setAttribute("refY","3");
    markerHead.setAttribute("orient","auto");
    markerHead.setAttribute("markerUnits","strokeWidth");

    let arrowPath = newSvgElmnt(markerHead,"path");
    arrowPath.setAttribute("d","M0,0 L0,6 L9,3 z");
    arrowPath.setAttribute("fill",color);

    line.setAttribute("marker-end","url(#arrowhead)");


}


function drawPolygon(SVG,polygonId, points, color = "black"){
  // if the polygon already exists, we'll update it
  // if it exists and the points given doesnt match the number of points in the polygon, we'll throw an error
  let polygon = SVG.querySelector(`#${polygonId}`);
  if (polygon){
      if (polygon.points.length != points.length){
          throw new Error(`Invalid number of points for polygon ${polygonId}`);
      }
      let pointsStr = points.map(p => p.join(",")).join(" ");
      polygon.setAttribute("points", pointsStr);
  }
  else{
      let polygon = newSvgElmnt(SVG,"polygon");
      let pointsStr = points.map(p => p.join(",")).join(" ");
      polygon.setAttribute("points", pointsStr);
      polygon.setAttribute("fill", color);
      polygon.setAttribute("id", polygonId);
  }

  return polygon;
}

function drawParticles(svgd3, xs, radii, groupParticlesId,color = "red", debug = false) {
    // Create or select the group for all particles
    // if radii is a number, we'll make it an array of the same length as xs with the same value
    if (typeof radii === "number"){
        radii = Array(xs.length).fill(radii);
    }
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

    if (debug){
    particles.select("text")
        .attr("x", (d, i) => d.x)
        .attr("y", (d, i) => d.y)
        .text((d, i) => `${i}`);
    }
    // Remove old particles (exit)
    particles.exit().remove();
  }




   //lets try to write drawDistanceConstraints without the annoying d3, lets use plain js and svg
   //we'll keep the signature for compatibility
    function drawDistanceConstraints(svg, xarray, constraints_idxs,constraints_distances,
                                  groupId, color = "black",labels_segments = [], debug = false){


      //Lets create a group where we'll place all elements, for easy removal

      let group = svg.querySelector(`#${groupId}`);
      if (group == null){
          group = newSvgElmnt(svg,"g");
          group.setAttribute("id",groupId);
      }

      //if some element exist we update the properties, else we create them
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
          let line = group.querySelector(`#line_${groupId}_${i}`);
          if (line == null){
              line = newSvgElmnt(group,"line");
              line.setAttribute("id",`line_${groupId}_${i}`);
              }
          line.setAttribute("x1", x1);
          line.setAttribute("y1", y1);
          line.setAttribute("x2", x2_imposed);
          line.setAttribute("y2", y2_imposed);
          line.setAttribute("stroke", color);
          //lets add the label if it exists

          if (debug){
          if (labels_segments.length > 0){
              if (labels_segments.length != constraints_idxs.length){
                  throw new Error("labels_segments must have the same length as constraints_idxs");
              }
              let xmid = (x1 + x2)/2;
              let ymid = (y1 + y2)/2;
              let label = labels_segments[i];
              let text = group.querySelector(`#label_${groupId}_${i}`);
              if (text == null){

                  text = newSvgElmnt(group,"text");
                  text.setAttribute("id",`label_${groupId}_${i}`);
              }
              text.setAttribute("x", xmid);
              text.setAttribute("y", ymid);
              text.textContent = label;
              text.setAttribute("fill", color);
              text.setAttribute("font-size", "18px");
              text.setAttribute("text-anchor", "middle");
          }
        }
      }

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


//again, lets write drawPinConstraints without d3, lets use plain js and svg

function drawPinConstraints(svg, xarray, constraints_idxs,
                            pinpoints ,distances, groupId, color = "black",labels_segments = [], debug = false){
  //very similar to drawDistanceConstraints but now there is only 1 particle involved
  // xarray [Nparticles,2] array, constraints_idxs [Nconstraints] pinpoints is a [Nconstraints,2]
  //the other coordinates correspond to the pin point
  // Lets throw an error if different lengths of pinpoints and constraints_idxs

  if (constraints_idxs.length != pinpoints.length){
    throw new Error("constraints_idxs and pinpoints must have the same length");
  }

  // Lets create a group where we'll place all elements, for easy removal

  let group = svg.querySelector(`#${groupId}`);
  if (group == null){

      group = newSvgElmnt(svg,"g");
      group.setAttribute("id",groupId);

  }

  //if some element exist we update the properties, else we create them

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
      let line = group.querySelector(`#line_${groupId}_${i}`);
      if (line == null){
          line = newSvgElmnt(group,"line");
          line.setAttribute("id",`line_${groupId}_${i}`);
      }
      line.setAttribute("x1", x1);
      line.setAttribute("y1", y1);
      line.setAttribute("x2", x2_imposed);
      line.setAttribute("y2", y2_imposed);
      line.setAttribute("stroke", color);
      //lets add the label if it exists

      if (debug){
        if (labels_segments.length > 0){

            if (labels_segments.length != constraints_idxs.length){

                throw new Error("labels_segments must have the same length as constraints_idxs");

            }

            let xmid = (x1 + x2)/2;
            let ymid = (y1 + y2)/2;
            let label = labels_segments[i];
            let text = group.querySelector(`#label_${groupId}_${i}`);
            if (text == null){

                text = newSvgElmnt(group,"text");
                text.setAttribute("id",`label_${groupId}_${i}`);

            }

            text.setAttribute("x", xmid);
            text.setAttribute("y", ymid);

            text.textContent = label;
            text.setAttribute("fill", color);
            text.setAttribute("font-size", "18px");
            text.setAttribute("text-anchor", "middle");
        }

      }

  }

  // additionally lets draw pin point with a cross, its coordinates are x2,y2
  // using the function drawCross, we'll pass svgd3 instead of svg for now
  for (let i = 0; i<constraints_idxs.length;i++){
    let x2 = pinpoints[i][0];
    let y2 = pinpoints[i][1];
    let id_cross = `pin_cross_${groupId}_${i}`;
    let svgd3 = d3.select(svg);
    drawCross(svgd3, x2, y2, 10, id_cross, color = "black");
  }

}





  function drawArrowd3(svgd3, x1, y1, x2, y2, arrowId, class_name = "arrow", color = "black") {
    // Check if the arrow group already exists
    let arrow = svgd3.select(`#${arrowId}`);
    if (arrow.empty()) {
        // Arrow doesn't exist, create new
        arrow = svgd3.append("g").attr("id", arrowId);
        // Create line
        arrow.append("line")
            .attr("id", `${arrowId}_line`);
        // Create defs and marker
        let defs = arrow.append("defs").append("marker")
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
            .attr("fill", color);
    }
    // Update line attributes
    arrow.select(`#${arrowId}_line`)
        .attr("x1", x1)
        .attr("y1", y1)
        .attr("x2", x2)
        .attr("y2", y2)
        .attr("stroke-width", 2)
        .attr("stroke", color)
        .attr("class", class_name)
        .attr("marker-end", `url(#${arrowId}_head)`);

    // Update marker color if needed
    arrow.select(`#${arrowId}_head path`)
        .attr("fill", color);
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
          console.debug("Plotting for the first time, creating axis for div_id", this.div_id)
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
          console.debug("Color not provided, using color", color)
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
  


function drawForces(STATE, svg, scaleArrows = 1.0){

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
      let ftotal = STATE.total_forces;
      let factuators = STATE.actuator_forces;
      let xs = STATE.xs;
      let sc = scaleArrows;
      for (let i = 0; i < xs.length; i++){
          let [x,y] = xs[i];
          let [fx, fy] = fext[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_fext_" + i, "arrow force fext", "blue");
          [fx, fy] = fconstraint[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_fconstraint_" + i, "arrow force fconstraint", "green");
          [fx, fy] = fcontact[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_fcontact_" + i, "arrow force fcontact", "red");
          [fx, fy] = ffriction[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_ffriction_" + i, "arrow force ffriction", "purple");
          [fx, fy] = ftotal[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_ftotal_" + i, "arrow force ftotal", "black");
          [fx, fy] = factuators[i];
          drawArrowd3(svg, x, y, x+sc*fx, y+sc*fy, "arrow_factuators_" + i, "arrow force factuators", "orange");
      }
  }


  function drawConstraints(STATE, svg, constraintsGroupId = "constraints", debug = false){
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
      //lets create labels, with the segment indices is basically an 0--nDistanceConstraints-1 array
      let labels_segments = Array.from(Array(constraints_idxs1.length).keys()).map(i => `Segd${i}`);
      let constraintsDistId = constraintsGroupId + "_distance";
      drawDistanceConstraints(svg, xs, constraints_idxs1, distances1, constraintsDistId, "black", labels_segments,debug);
    }
    if (!isConsPinEmpty){
      let constraints_pin = STATE.constraints_pin;
      let constraints_idxs2 = constraints_pin.map(constraint => constraint[0]);
      let pinpoints = constraints_pin.map(constraint => constraint[1]);
      let distances2 = constraints_pin.map(constraint => constraint[2]);
      //now the labels go form nDistanceConstraints to nDistanceConstraints + nPinConstraints - 1
      let labels_segments = Array.from(Array(constraints_idxs2.length).keys()).map(i => `Segp${i}`);
      let constraintsPinId = constraintsGroupId + "_pin";
      drawPinConstraints(svg, xs, constraints_idxs2, pinpoints, distances2, constraintsPinId, "black", labels_segments,debug);
    }
    
}



class PolygonDrawer{
    /*
    * Here we pack all the functions to draw polygons in the svg, to take care of normals, transformations etc
    * It is intended to be the visual representation of the Polygon object
    */
    constructor(svg,VerticesCanvas,polId,polygonsGroupId = "polygonsGroup", color = "blue", debug = false){
      this.polygonObj = new Polygon(VerticesCanvas, polId);
      this.polygonsGroupId = polygonsGroupId;
      //polygons in the svg will have the id ${polygonGroupId}_polygon_${i} where i is the index/id of the polygon
      // any property of the polygon will be ${polygonGroupId}_polygon_${i}_attribute 
      //for instance ${polygonGroupId}_polygon_${i}_normal_${j} where j is the index of the normal
      //lets draw the polygon
      this.debug = debug;
      this.draw(svg, color);
      

    }

    draw(svg, color = "blue", scale_normal = 20.0){

        //lets retrieve/create a polygon group that will contain the polygon and the normals and any other property
        //lets use vanilla js 
        let polygonGroupInnerId = this.polygonsGroupId + "_polygon_" + this.polygonObj.id;
        let polygonGroup = document.getElementById(polygonGroupInnerId);
        if (polygonGroup === null){
          polygonGroup = newSvgElmnt(svg,"g")
          polygonGroup.id = polygonGroupInnerId;
        }

        //lets draw the polygon
        let points = this.polygonObj.vertices;
        drawPolygon(polygonGroup,polygonGroup.id+"_polygon", points, "blue")

        //Lets display the id of the polygon at the center of the polygon, plain js

        if (this.debug){
            let centroid = math.mean(points,0)
            let textId = polygonGroup.id + "_id";
            let text = polygonGroup.querySelector(`#${textId}`);
            if (text === null){
              text = newSvgElmnt(polygonGroup,"text");
              text.id = textId;
              text.setAttribute("x",centroid[0]);
              text.setAttribute("y",centroid[1]);
              text.setAttribute("fill","black");
              text.setAttribute("font-size","25px");
              text.setAttribute("text-anchor","middle");
              text.setAttribute("alignment-baseline","middle");
              text.textContent = this.polygonObj.id;
            }

            
            //lets draw the normals
            this.drawNormals(polygonGroup, "black", scale_normal);

      }

    }

    drawNormals(svg, color = "black",scale_normal = 1.0){
      //lets draw the normals of the polygon
      
      let normals = this.polygonObj.getNormals();
      let polygonGroupInnerId = this.polygonsGroupId + "_polygon_" + this.polygonObj.id;
      let polygonGroup = document.getElementById(polygonGroupInnerId);
      if (polygonGroup === null){
        throw new Error("Polygon group must be created before drawing normals");
      }
      for (let i = 0; i < normals.length; i++){
        let normal = normals[i];
        let normalId = polygonGroupInnerId + "_normal_" + i;
        let midPoint = this.polygonObj.getEdgeMidpoint(i);
        //lets remove the normal if it already exists
        let previousNormal = document.getElementById(normalId);
        

        drawArrowd3(d3.select(polygonGroup), midPoint[0],midPoint[1],
        midPoint[0] - normal[0]*scale_normal, midPoint[1] - normal[1]*scale_normal, normalId, "normal", color);
        //The minuss sign is because the 90º rotation in canvas coordinates (y grows down) 
        //Lets add a text with the normal index
        // we check if the text already exists and if it does we update it
        let textId = normalId + "_text";
        let text = polygonGroup.querySelector(`#${textId}`);
        if (text === null){
          text = newSvgElmnt(polygonGroup,"text");
          text.id = textId;
          text.setAttribute("x",midPoint[0] - normal[0]*scale_normal);
          text.setAttribute("y",midPoint[1] - normal[1]*scale_normal);
          text.setAttribute("fill","black");
          text.setAttribute("font-size","12px");
          text.setAttribute("text-anchor","middle");
          text.setAttribute("alignment-baseline","middle");
          text.textContent = i;
        }
        else{
          text.setAttribute("x",midPoint[0] - normal[0]*scale_normal);
          text.setAttribute("y",midPoint[1] - normal[1]*scale_normal);
        }


      }

    

    }
  }
class CollisionDrawer{
  /*
  * Class to draw the collisions in the svg
  *
  */
        constructor(svg, particles_group_id, polygonsGroupId){
            this.svg = svg;
            this.particles_group_id = particles_group_id;
            this.polygons_group_id = polygonsGroupId;
            //particles id are expected to be ${particles_group_id}_particle_${i} where i is the index of the particle
            //polygons id are expected to be ${polygons_group_id}_polygon_${i} where i is the index of the polygon
            //any property of the polygon is expected to be:
            //${polygons_group_id}_polygon_${idx}_attribute for instance ${polygons_group_id}_polygon_${idx}_normal

            
        }
        

        drawCollisions(new_collisions, resolved_collisions){
            //new_collisions and resolved_collisions are arrays of type [particle_index,edge_index,polygon_idx]
            //      

            for ( let collision of new_collisions){
                this._enterCollision(collision);
            }

            for (let collision of resolved_collisions){
                this._exitCollision(collision);
            }

        }

        _enterCollision(collision){
            //We'll change the edge normal to red and the particle to red when in contact
            let [particle_index,edge_index,polygonIdx] = collision;
            let particle = this._getParticleByIdx(particle_index);
            let normal = this._getNormalByIdxPolygon(polygonIdx, edge_index);

            particle.attr("fill","red");
            normal.attr("stroke","red");

            

        }

        _exitCollision(collision){
            //We'll change the edge normal to black and the particle to green when in contact
            let [particle_index,edge_index,polygonIdx] = collision;
            let particle = this._getParticleByIdx(particle_index);
            let normal = this._getNormalByIdxPolygon(polygonIdx, edge_index);
            particle.attr("fill","green");
            normal.attr("stroke","black");
        }

      
      _getParticleByIdx(idx){
        let particle = d3.select(`#${this.particles_group_id}_particle_${idx}`);
        //lets throw an error if the particle is not found
        if (particle.empty()){
          let query = `#${this.particles_group_id}_particle_`;
          throw new Error(`Particle with index ${idx} not found, query: ${query}`);
        }
        return particle;
      }
      _getNormalByIdxPolygon(polIdx, edgeIdx){
        let normal = d3.select(`#${this.polygons_group_id}_polygon_${polIdx}_normal_${edgeIdx}`);
        //lets throw an error if the polygon normal is not found
        
        if (normal.empty()){
          throw new Error(`Normal with index ${edgeIdx} not found in polygon ${polIdx}`);
        }
        return normal;

      }

    }

  function drawSpring(svg,point1, point2,
        restLength,springGroupId, nSegments=10) {
    const dx = point2[0] - point1[0];
    const dy = point2[1] - point1[1];
    const dist = Math.sqrt(dx * dx + dy * dy);
    const dir = [dx / dist, dy / dist]; // Unit direction vector
    const perp = [-dir[1], dir[0]]; // Perpendicular vector
    const width = 0.1 * restLength; // Spring width

    //lets remove the previous spring if it exists
    let previousSpring = document.getElementById(springGroupId);
    if (previousSpring){
    svg.removeChild(previousSpring);
    }
    //lets create a group for the spring
    let springGroup = newSvgElmnt(svg,"g");
    springGroup.setAttribute("id",springGroupId);



    for (let i = 0; i < nSegments; i++) {
    const fraction = i / nSegments;
    const nextFraction = (i + 1) / nSegments;
    const offsetDirection = i % 2 == 0 ? 1 : -1; // Alternate direction

    // Adjust start and end points for the first and last segments
    const startPoint = i === 0 ? [point1[0], point1[1]] : [
    point1[0] + fraction * dx + offsetDirection * perp[0] * width,
    point1[1] + fraction * dy + offsetDirection * perp[1] * width,
    ];
    const endPoint = i === nSegments - 1 ? [point2[0], point2[1]] : [
    point1[0] + nextFraction * dx - offsetDirection * perp[0] * width,
    point1[1] + nextFraction * dy - offsetDirection * perp[1] * width,
    ];

    // Draw this segment
    let line = newSvgElmnt(springGroup, "line");
    line.setAttribute("x1", startPoint[0]);
    line.setAttribute("y1", startPoint[1]);
    line.setAttribute("x2", endPoint[0]);
    line.setAttribute("y2", endPoint[1]);
    line.setAttribute("stroke", "black");
    line.setAttribute("stroke-width", "2");
    }


}



function particleInteractionMouse(SVG, particle_id, canvas2Model, mouseMoveCallbacks = [], mouseUpCallbacks = []){
  //this function is intended to execute the provided callbacks when the particle is clicked
  // all mouseMoveCalbacks must have the signature callback(mouseCoordsModel, particle_idx)
  // all mouseUpCallbacks must have the signature callback()
  // with mouseup we remove the mousemove and mouseup event listeners
  let particle = document.getElementById(particle_id);
  //lets throw an error if the particle is not found
  if (particle === null){
      throw new Error("Particle not found with id " + particle_id);
  }
  //lets retrieve the particle_idx ( it must be the last string after _ in the id)
  let split_id = particle_id.split("_");
  //parseInt
  let particle_idx = parseInt(split_id[split_id.length-1]);
  //if not a number, throw an error
  if (isNaN(particle_idx)){
      throw new Error("Particle id must end with a number");
  }
  particle.addEventListener("mousedown", function(event){
      event.preventDefault();

      function mouseMoveCallback(event){
          

          let [xabs, yabs] = [event.clientX, event.clientY];
          SVG.lastClientX = xabs;
          SVG.lastClientY = yabs; //super sloppy fix to manually update mousemove event
          let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);
          let [xnewCanvas, ynewCanvas] =[xsvg, ysvg]
          let [xnewModel, ynewModel] = canvas2Model([[xnewCanvas, ynewCanvas]])[0];
          let mouseCoordsModel = [xnewModel, ynewModel];
          for (let callback of mouseMoveCallbacks){
              callback(mouseCoordsModel, particle_idx);
          }
      }
      function mouseUpCallback(event){
          SVG.removeEventListener("mousemove", mouseMoveCallback);
          SVG.removeEventListener("mouseup", mouseUpCallback);

          for (let callback of mouseUpCallbacks){
              callback();
          }
      }

      SVG.addEventListener("mousemove", mouseMoveCallback);
      SVG.addEventListener("mouseup", mouseUpCallback);

  })

}
     



function _mouseUpSpringCallback(drawer,STATE){
  STATE.springs2points = STATE.springs2points.filter((spring) => {
      return spring.length != 5 || spring[4] != "mouseSpring";
  });

  drawer.removeSpringMouse();
}


function _mouseMoveSpringCallback(drawer,STATE,mouseSpringRestLength,
  mouseSpringStiffness,
  mouseCoords, particle_idx){
  //we attach a mouseSpring text tag to the last entry of springs2points so we can identify it for update and removal
  // we iter over the springs2points, looking for an array with 5 elements, the last one being a mouseSpring tag
  // if we find it we update the spring to the mouse position
  // if we don't find it we create a new spring
  //




  let spring = STATE.springs2points.find((spring) => {
      return spring.length == 5 && spring[4] == "mouseSpring";
  });
  if (spring){
      spring[1] = mouseCoords;
  }else{
      spring = [particle_idx,mouseCoords,mouseSpringStiffness,
        mouseSpringRestLength,"mouseSpring"];
      STATE.springs2points.push(spring);
  }

  let [pIdx,mouseCoordsModel,_1,restLengthModel,_2] = spring;
  drawer.drawSpringMouse(STATE, pIdx,mouseCoordsModel,restLengthModel);




}



class Drawer{
      /*
      * Lets make a recopilation of all the state properties drawn by the Drawer
      * xs: [Nparticles,2] array with the positions of the particles
      * vs: [Nparticles,2] array with the velocities of the particles
      * external_forces: [Nparticles,2] array with the external forces acting on the particles
      * constraint_forces: [Nparticles,2] array with the constraint forces acting on the particles
      * contact_forces: [Nparticles,2] array with the contact forces acting on the particles
      * friction_forces: [Nparticles,2] array with the friction forces acting on the particles
      * constraints_distance: [Nconstraints,3] array with the distance constraints, [idx1,idx2,distance]
      * constraints_pin: [Nconstraints,3] array with the pin constraints, [idx,[x,y],distance]
      * new_collisions: [Ncollisions,3] array with the new collisions, [particle_index,edge_index,polygon_idx]
      * resolved_collisions: [Ncollisions,3] array with the resolved collisions, [particle_index,edge_index,polygon_idx]
      * polygons: [[verticesPol1],[verticesPol2],...] array with the vertices of the polygons
      * springs: [[idx1,idx2,springConstant,restLength],...] array with the springs
      * springs2points: [[idx,point,springConstant,restlength]]
      * 
      */
      constructor(svg,options){

        let defaultOptions = {
          x_model_domain: [-1,1],
          y_model_domain: [-1,1],
          origin_model: [0,0],
          preffixId: "",
          debug: false
        };
        for (let key in options){
          defaultOptions[key] = options[key];
        }

        this.svg = svg;// html dom, not d3
        this.x_model_domain = defaultOptions.x_model_domain;
        this.y_model_domain = defaultOptions.y_model_domain;
        this.origin_model = defaultOptions.origin_model;
        this.preffixId = defaultOptions.preffixId? preffixId : svg.id;
        this._firstMakeInteractiveCall = true;
        this.debug = defaultOptions.debug;
        
        //if prefixId is provided we'll use it, if not we get the id of the svg
      }

      _checkInitState(STATE){
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

        return STATE;
      }

      getCanvasDimensions(){
        let width = this.svg.clientWidth;
        let height = this.svg.clientHeight;
        //lets throw an error if width or height are 0
        if (width === 0 || height === 0){
          throw new Error("Width or height of the svg is 0");
        }

        let origin_canvas = [width/2,height/2];// canvas point where the origin will be placed
        let scale_factor_x = width/(this.x_model_domain[1] - this.x_model_domain[0]);
        let scale_factor_y = -height/(this.y_model_domain[1] - this.y_model_domain[0]);//negative because in svg y grows down

        return [width,height,origin_canvas,scale_factor_x,scale_factor_y];
      }
      model2Canvas(array, shift = false){
        //lets check type, array is [N,2] model data is expected to be O(1)
        //shift is an optional parameter because coord points will be shifted by the origin but forces will not ( second derivative of the position)
        //lets turn the array into array if it is a mathjs matrix

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

      _model2CanvasSpring(STATE,spring){
        let [idx,pinpoint,springConstant,restLength] = spring;
        let point1 = STATE.xs[idx];
        let dx = pinpoint[0] - point1[0];
        let dy = pinpoint[1] - point1[1];
        let current_distance = Math.sqrt(dx**2 + dy**2);
        let [point1Canvas,pinpointCanvas] = this.model2Canvas([point1,pinpoint],true);
        let distanceCanvas = this._calculateDistance(point1Canvas,pinpointCanvas);
        let ratio = distanceCanvas/current_distance;
        let restLengthCanvas = ratio * restLength;

        return [idx,pinpointCanvas,springConstant,restLengthCanvas];
      }

      state2Canvas(STATE){
        // Here we transform the different magnitudes in STATE to canvas coordinates, we must generate a twin canvasSTATE object
        //lets start by deep copying the STATE object
        // positions are transformed with  shift = true, forces and velocities are transformed with shift = false

        let canvasSTATE = JSON.parse(JSON.stringify(STATE));
        //lets turn everything into arrays if they are not
        
        canvasSTATE.xs = this.model2Canvas(STATE.xs,true);
        canvasSTATE.vs = this.model2Canvas(STATE.vs,false);
        canvasSTATE.external_forces = this.model2Canvas(STATE.external_forces,false);
        canvasSTATE.constraint_forces = this.model2Canvas(STATE.constraint_forces,false);
        canvasSTATE.contact_forces = this.model2Canvas(STATE.contact_forces,false);
        canvasSTATE.friction_forces = this.model2Canvas(STATE.friction_forces,false);
        canvasSTATE.total_forces = this.model2Canvas(STATE.total_forces,false);
        canvasSTATE.actuator_forces = this.model2Canvas(STATE.actuator_forces,false);

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


        //lets transform the polygons
        let polygonsVertices = STATE.polygons;
        let canvasPolygonsVertices = polygonsVertices.map(vertices => this.model2Canvas(vertices,true));

        canvasSTATE.polygons = canvasPolygonsVertices;


        // we have to transform the spring points and restLengths
        //spring forces come in 2 types, between particles and between particle and point, springs: [[idx1,idx2,springConstant,restLength],...]
        //springs2points: [[idx1,[x,y],springConstant,restLength],...]
        //to transform restLengths we get the ratio of canvasDistance/modelDistance for the \vec{dr} and multiply it by the restLength
        // first springs
        if (STATE.hasOwnProperty("springs")){
          let canvasSprings = STATE.springs.map(spring => {
            let [idx1,idx2,springConstant,restLength] = spring;
            let point1 = STATE.xs[idx1];
            let point2 = STATE.xs[idx2];
            let dx = point2[0] - point1[0];
            let dy = point2[1] - point1[1];
            let current_distance = Math.sqrt(dx**2 + dy**2);
            let [point1Canvas,point2Canvas] = this.model2Canvas([point1,point2],true);
            let distanceCanvas = this._calculateDistance(point1Canvas,point2Canvas);
            let ratio = distanceCanvas/current_distance;
            let restLengthCanvas = ratio * restLength;
            
            return [idx1,idx2,springConstant,restLengthCanvas];
          });

          canvasSTATE.springs = canvasSprings;
        }

        if (STATE.hasOwnProperty("springs2points")){

          /*
          let canvasSprings2Points = STATE.springs2points.map(spring => {
          //we'll skip the mouseSpring, its drawing is managed by the mouseMoveSpringCallback
          

            let [idx,pinpoint,springConstant,restLength] = spring;
            let point1 = STATE.xs[idx];
            let dx = pinpoint[0] - point1[0];
            let dy = pinpoint[1] - point1[1];
            let current_distance = Math.sqrt(dx**2 + dy**2);
            let [point1Canvas,pinpointCanvas] = this.model2Canvas([point1,pinpoint],true);
            let distanceCanvas = this._calculateDistance(point1Canvas,pinpointCanvas);
            let ratio = distanceCanvas/current_distance;
            let restLengthCanvas = ratio * restLength;

            
            return [idx,pinpointCanvas,springConstant,restLengthCanvas];
          });
          */

          //lets filter to avoid the mouseSpring
          let canvasSprings2Points = STATE.springs2points.filter(spring => {
            return spring.length != 5 || spring[4] != "mouseSpring";
          }).map(spring => {

              
              return this._model2CanvasSpring(STATE,spring);
            })

          canvasSTATE.springs2points = canvasSprings2Points;
        }






        return canvasSTATE;


      }


      drawPolygons(svg, polygonsVertices,
                         polygonsGroupId, color = "blue"){
        for (let [index,vertices] of polygonsVertices.entries()){
          new PolygonDrawer(svg,vertices,index,polygonsGroupId,color,this.debug);
      }
      }

      drawCollisions(svg, new_collisions,
                     resolved_collisions,
                     particlesGroupId, polygonsGroupId){
        let collisionDrawer = new CollisionDrawer(svg, particlesGroupId, polygonsGroupId);
        collisionDrawer.drawCollisions(new_collisions, resolved_collisions);
      }

      drawSprings(svg, points1, points2, restLengths, springsGroupId){
        //if points1 is empty we skip the drawing
        if (points1.length !== 0){
        // lets throw an error if the points1, points2 and restLengths are not of the same length
        if (points1.length !== points2.length || points1.length !== restLengths.length){
          throw new Error("points1, points2 and restLengths must be of the same length");
        }
        //lets create the group of springs
        let springsGroup = document.getElementById(springsGroupId);
        if (springsGroup){
          svg.removeChild(springsGroup);
        }
        springsGroup = newSvgElmnt(svg,"g");
        springsGroup.id = springsGroupId;
        for (let i = 0; i < points1.length; i++){
          let springId = springsGroupId + "_spring_" + i;
          drawSpring(springsGroup,points1[i],points2[i],restLengths[i],springId);
        }
      }


      }

      drawSpringMouse(STATE, pIdx,mouseCoordsModel,restLengthModel){


        let springId = `spring_mouse`;
        let [idx,pinpointCanvas,_,restLengthCanvas] = this._model2CanvasSpring(STATE,[pIdx,mouseCoordsModel,-1,restLengthModel]);
        let point1Canvas = this.model2Canvas([STATE.xs[pIdx]],true)[0];
        drawSpring(this.svg,point1Canvas,pinpointCanvas,restLengthCanvas,springId);

      }

      removeSpringMouse(id = "spring_mouse"){
        let spring = document.getElementById(id);
        if (spring){
          this.svg.removeChild(spring);
        }

      }


      drawState(STATE, particlesGroupId = "particlesGroup",
                polygonsGroupId = "polygonsGroup",
                 VISUAL_OPTIONS = {}, springsGroupId = "springsGroup", constraintsGroupId = "constraintsGroup"){

        // lets get the id of this.svg to use it as prefix for everything drawn here
        //to get the id of an svg element we use the id property
        let preId = this.preffixId;
        let particlesGroupId_full = preId + "_" + particlesGroupId;
        let polygonsGroupId_full = preId + "_" + polygonsGroupId;
        let springsGroupId_full = preId + "_" + springsGroupId;
        let constraintsGroupId_full = preId + "_" + constraintsGroupId;

        this.particlesGroupId = particlesGroupId_full;
        this.polygonsGroupId = polygonsGroupId_full;
        this.springsGroupId = springsGroupId_full;
        this.constraintsGroupId = constraintsGroupId_full;


        



        let defaultVisualOptions = {
          particle_density: 0.0005,
          scale_arrows_forces: 1.0,
          color_particles: "red",
          color_polygons: "blue",
        }
        //we iter over defaultVisualOptions, and if the key is not in VISUAL_OPTIONS we add it
        for (let key in defaultVisualOptions){
          if (!VISUAL_OPTIONS.hasOwnProperty(key)){
            VISUAL_OPTIONS[key] = defaultVisualOptions[key];
          }
        }

        
        //lets check the state and initialize missing properties
        STATE = this._checkInitState(STATE);
        let s = STATE;
        let svgd3 = d3.select(this.svg);
        let particle_density = VISUAL_OPTIONS["particle_density"];
        let radii = s.masses.map(mass => Math.sqrt(mass/particle_density));
        let scaleArrowsForces = VISUAL_OPTIONS["scale_arrows_forces"];
        let canvasSTATE = this.state2Canvas(s);

        //lets draw particles
        drawParticles(svgd3, canvasSTATE.xs, radii, particlesGroupId_full, VISUAL_OPTIONS["color_particles"], this.debug);
        this.drawPolygons(this.svg, canvasSTATE.polygons, polygonsGroupId_full, VISUAL_OPTIONS["color_polygons"]);
        drawConstraints(canvasSTATE, this.svg, constraintsGroupId_full, this.debug);
        //lets draw forces
        if(this.debug){
        drawForces(canvasSTATE, svgd3, scaleArrowsForces);
        //lets draw constraints
        


        //lets draw the collisions
        this.drawCollisions(this.svg, canvasSTATE.new_collisions,
                            canvasSTATE.resolved_collisions,
                              particlesGroupId_full, polygonsGroupId_full);
        }
        //lets draw the springs
        //we first extract the points and restLengths from both springs and springs2points, together with the restlengths, and we concatenate both sets
        let points1 = canvasSTATE.springs.map(spring => canvasSTATE.xs[spring[0]]);
        let points2 = canvasSTATE.springs.map(spring => canvasSTATE.xs[spring[1]]);

        let restLengths = canvasSTATE.springs.map(spring => spring[3]);
        let points1_2 = canvasSTATE.springs2points.map(spring => canvasSTATE.xs[spring[0]]);
        let points2_2 = canvasSTATE.springs2points.map(spring => spring[1]);
        let restLengths_2 = canvasSTATE.springs2points.map(spring => spring[3]);
        points1 = points1.concat(points1_2);
        points2 = points2.concat(points2_2);
        restLengths = restLengths.concat(restLengths_2);
        this.drawSprings(this.svg, points1, points2, restLengths, springsGroupId_full);

        this.makeInteractive(STATE);
    

      }

      getParticlesGroup(){
        return document.getElementById(this.particlesGroupId);
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


      makeInteractive(STATE){

        /////DISGUSTING PATCH TO MANUALLY DISPATCH MOUSEMOVE EVENT AND UPDATE
        //MOUSE SPRING ATTACHMENT POINT
        
        let [lastClientX, lastClientY] = [this.svg.lastClientX, this.svg.lastClientY];

        if (typeof(lastClientX) === "number" && typeof(lastClientY) === "number"){
          let mousemoveEvent = new MouseEvent("mousemove", {
            view: window,
            bubbles: true,
            cancelable: true,
            clientX: lastClientX,
            clientY: lastClientY
          });
          this.svg.dispatchEvent(mousemoveEvent);
        }




        ////

        if (this._firstMakeInteractiveCall){
          this._firstMakeInteractiveCall = false;
        }
        else{
          return;
        }
        console.info(` Enable interactive particles `)
        let [mouseSpringRestLength, mouseSpringStiffness] = [STATE.mouseSpringRestLength, STATE.mouseSpringStiffness];

      //if null we default them to 10 and 0.1
      if (typeof(mouseSpringStiffness) !== "number"){
          mouseSpringStiffness = 5;

          console.warn("Mouse spring stiffness not found in STATE, defaulting");
      }

      if (typeof(mouseSpringRestLength) !== "number"){

          mouseSpringRestLength = 0.5;
          console.warn("Mouse spring rest length not found in STATE, defaulting");
      }


        let canvas2ModelCallback = (pointCanvas) => {
          return this.canvas2Model(pointCanvas,true);

          }

        let mouseMoveSpringCallback = (mouseCoordsModel, particle_idx) => {
          _mouseMoveSpringCallback(this,STATE, mouseSpringRestLength,
            mouseSpringStiffness,
            mouseCoordsModel, particle_idx);
        }

        let mouseUpSpringCallback = () => {
          _mouseUpSpringCallback(this,STATE);
        }

        let particlesGroup = this.getParticlesGroup();
        if (particlesGroup){
          let particles = particlesGroup.querySelectorAll(".particle");
          particles.forEach((particle,idx) => {
              particleInteractionMouse(this.svg, particle.id,
                 canvas2ModelCallback, [mouseMoveSpringCallback], [mouseUpSpringCallback]);
          }
          );
        }
      }


 
      


      remove(){
        //lets remove all the groups in the svg
        let groups = [this.particlesGroupId, this.polygonsGroupId, this.springsGroupId, this.constraintsGroupId+"_distance", this.constraintsGroupId+"_pin"];

        for (let group of groups){
          let groupElmnt = document.getElementById(group);
          if (groupElmnt){
            this.svg.removeChild(groupElmnt);
          }
        }
      }
    

    }
export {drawCircles, drawLines,  drawArrow,drawPolygon,
   drawParticles, drawDistanceConstraints,
   drawPinConstraints, drawCross, drawForces,
   drawConstraints,drawSpring,
   particleInteractionMouse,
    CollisionDrawer,
    Plot, Drawer,PolygonDrawer};