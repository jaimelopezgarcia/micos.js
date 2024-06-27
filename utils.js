

/**
 * Slices a given matrix according to the provided slices.
 *
 * @param {Matrix} matrix - The matrix to slice.
 * @param {Array} slices - An array of slices. Each slice can be a single integer (for a specific index) or an array of two integers (for a range of indices).
 * @returns {Matrix} The sliced matrix.
 *
 * @example
 * // Create a 3x3 matrix
 * let m = math.matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
 *
 * // Slice the matrix to get the first two rows and the first column
 * let sliced = sliceMatrix(m, [[0, 2], 0]);
 * // sliced is now a 2x1 matrix: [[1], [4]]
 */
function sliceMatrix(matrix, slices) {
    let shape = matrix.size();
    let ranges = [];
  
    for (let i = 0; i < shape.length; i++) {
      let slice = slices[i];
      if (Array.isArray(slice)) {
        let start = slice[0];
        let end = slice[1];
        if (end == -1) {
          end = shape[i];
        }
        ranges.push(math.range(start, end));
      } else {
        ranges.push(slice); // is an int
      }
    }
  
    let sliced = math.subset(matrix, math.index(...ranges));
    return sliced;
  }

  function prettyPrintState(state) {
    var output = "";
  
    var list_keys = ["xs", "vs", "external_forces",
                    "constraint_forces", "J",
                    "lagrange_multipliers", "dot_J",
                    "constraint_vals"];
    for (var key of list_keys) {
  
      if (state.hasOwnProperty(key)) {
        output += key + ":\n";
        var value = state[key];
        if (Array.isArray(value) && Array.isArray(value[0])) { // if value is a 2D array
          for (var i = 0; i < value.length; i++) {
            
            output += "  [" + value[i].map(function(num){return num.toFixed(3)}).join(", ") + "]\n";
          }
        } else if (Array.isArray(value)) { // if value is a 1D array
          output += "  [" + value.map(function(num){return num.toFixed(3)}).join(", ") + "]\n";
        } else { // if value is not an array
          output += "  " + value + "\n";
        }
      }
    }
    return output;
  }


  function displayStateText(state,div_id){
    let text = prettyPrintState(state);
    text = text.replace(/\n/g, "<br/>");
    let debug_div = document.getElementById(div_id);
    debug_div.innerHTML = text;
  }



function capDecimalsArray(array, decimals = 3) {

	
    if (!Array.isArray(array)) {
      throw new Error("Input is not an array");
    };
    let dimensions = getArrayDimensions(array);
    
    if ((dimensions.length!= 1) && (dimensions.length!=2)){
      throw new Error(`Only 1 or 2 dimensions not ${dimensions}`);
    };
    
  
    if (Array.isArray(array[0])) { // if array is 2D
      return array.map(function(subArray) {
        if (!Array.isArray(subArray)) {
          throw new Error("Input is not a 2D array");
        }
        return subArray.map(function(num) {
          return num.toFixed(decimals);
        });
      });
    } else { // if array is 1D
      return array.map(function(num) {
        return num.toFixed(decimals);
      });
    }
  }
  
  function getArrayDimensions(array) {
    let dimensions = [];
    while (Array.isArray(array)) {
      dimensions.push(array.length);
      array = array[0];
    }
    return dimensions;
  }


  function displayObjectInDiv(obj, div_output_id, id_prefix = ""){
    /*
    Display the object in a div with the id div_output_id
    Info displayed in a new paragraph in the parent div with div_output_id with id_prefix_key
    */
  
    //lets throw an error is obj is not an object or a Map
    if (typeof obj !== 'object' && !(obj instanceof Map)){
      throw new Error("Input is not an object or a Map");
    }
    let parentDiv = document.getElementById(div_output_id)
    
    // Use Object.entries for objects and Map.prototype.entries for Maps
    let entries = obj instanceof Map ? obj.entries() : Object.entries(obj);
  
    for (let [key, value] of entries){
      let idchild = `#${id_prefix}-${key}`
      let childElement = parentDiv.querySelector(idchild)
      if (childElement){
        let valueStr = JSON.stringify(value);
        childElement.textContent = `${key}: ${valueStr}`;
      }
      else{
        let newP = document.createElement("p");
        if (id_prefix !== ""){
          id_prefix = div_output_id;
        }
        newP.id = `${id_prefix}-${key}`;
        let valueStr = JSON.stringify(value);
        newP.textContent = `${key}: ${valueStr}`;
        parentDiv.appendChild(newP);
      }
    }
  }
  

function getSvgRelativeCoords(svg, xabs, yabs){
  let rect = svg.getBoundingClientRect();
  let [xoffset, yoffset] = [rect.x, rect.y];
  let [xsvg, ysvg] = [xabs - xoffset, yabs - yoffset];
  return [xsvg, ysvg];
}


function displayMouseCoords(svg,  canvas2ModelCallback = null,div_output_id = null) {
  svg.addEventListener("mousemove", function(event) {
    let [xabs, yabs] = [event.clientX, event.clientY];
    let [xsvg, ysvg] = getSvgRelativeCoords(svg, xabs, yabs);
    if (div_output_id) {
      let div = document.getElementById(div_output_id);
      div.textContent = `Mouse coordinates: (${xsvg.toFixed(2)}, ${ysvg.toFixed(2)})`;
    } else {
      // Check if the text element already exists
      let textElement = svg.querySelector("#mouse-coords");
      if (!textElement) {
        // If it doesn't exist, create it
        textElement = document.createElementNS("http://www.w3.org/2000/svg", "text");
        textElement.setAttribute("id", "mouse-coords");
        textElement.setAttribute("x", 10); // Adjust position as needed
        textElement.setAttribute("y", 20); // Adjust position as needed
        svg.appendChild(textElement);
      }
      // Update the text content
      if (canvas2ModelCallback) {
        let [xmodel, ymodel] = canvas2ModelCallback([ [xsvg, ysvg] ])[0];
        let strText = `Mouse coordinates: (${xsvg.toFixed(2)}, ${ysvg.toFixed(2)}) (${xmodel.toFixed(2)}, ${ymodel.toFixed(2)})`;
        textElement.textContent = strText;
      }
      else {
        textElement.textContent = `Mouse coordinates: (${xsvg.toFixed(2)}, ${ysvg.toFixed(2)})`;
      }
    }
  });
}


export { sliceMatrix, prettyPrintState, displayStateText, capDecimalsArray,
       getArrayDimensions, displayObjectInDiv, getSvgRelativeCoords, displayMouseCoords};