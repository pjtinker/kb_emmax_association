/* GWAS v1.0b */

var margin = { top: 50, right: 100, bottom: 50, left: 80 },
    outerWidth = 1050,
    outerHeight = 500,
    width = outerWidth - margin.left - margin.right,
    height = outerHeight - margin.top - margin.bottom;

var x = d3.scale.linear()
      .range([0, width]).nice().clamp(true);
    //.range([0, width]).nice();
    //.range([0, width]).clamp(true);

var y = d3.scale.linear()
    .range([height, 0]).nice();

var xCat = "CHR",
    yCat = "P",
    xLabel = "Chromosome";
    yLabel = "p Value";
    colorCat = "red";

var title;
var minBP = new Object();
var maxBP = new Object();

function load(url)
{
  $("#loading").css("visibility","visible");
  $("#file").val("");

  // remove any current plot
  d3.select("svg").remove();

  d3.tsv(url, function(error, data) {
    if (error) throw error;
  
    data.forEach(function(d) {
      d.CHR = +d.CHR;
  //  d.P = +d.P;
      d.P = - Math.log10(d.P);
      d.BP = +d.BP;

    // max/min of BP
    maxBP[d.CHR] = Math.max(maxBP[d.CHR] || 0, d.BP);
    minBP[d.CHR] = Math.min(minBP[d.CHR] || 0, d.BP);
    });

    setTimeout(function(){ $("#loading").css("visibility","hidden"); }, 1000);
  
      var xMax = d3.max(data, function(d) { return d[xCat]; }) * 1.05,
        xMin = d3.min(data, function(d) { return d[xCat]; }),
        xMin = xMin > 0 ? 0 : xMin,
        yMax = d3.max(data, function(d) { return d[yCat]; }) * 1.05,
        yMin = d3.min(data, function(d) { return d[yCat]; }),
        yMin = yMin > 0 ? 0 : yMin,
        nChrom = d3.max(data, function(d) { return d[xCat]; });

        
    var title = data.length + " samples in " + nChrom + " chromosomes";
  
    x.domain([xMin, xMax]);
    y.domain([yMin, yMax]);
  
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom")
        ///.tickFormat(d3.format("d"))
        .tickFormat(function (d) {if(d == 0) return ''; else return d3.format("d")(d); }) // hide axis origin!
        .outerTickSize(0)
        .tickSize(-height);

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left")
        .tickSize(-width);
  
    var color = d3.scale.category10();
  
    var tip = d3.tip()
        .attr("class", "d3-tip")
        .offset([-10, 0])
        .html(function(d) {
          return "SNP: " + d["SNP"] + "<br>" + xLabel + ": " + d[xCat] + "<br>" + "BP: " + d["BP"] + "<br>" + 
                 yCat + ": " + d[yCat];
        });
  
    var zoomBeh = d3.behavior.zoom()
        .x(x)
        .y(y)
        .scaleExtent([0, 500])
  //    .scaleExtent([0, 10])
        .on("zoom", zoom);
  
    var svg = d3.select("#scatter")
      .append("svg")
        .attr("width", outerWidth)
        .attr("height", outerHeight)
      .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .call(zoomBeh);

/*
    svg.selectAll(".tick")
    .filter(function (d) { return d === 0;  })
    .remove();
*/
  
    svg.call(tip);
  
    svg.append("rect")
        .attr("width", width)
        .attr("height", height);
  
    svg.append("g")
        .classed("x axis", true)
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis)
      .append("text")
        .classed("label", true)
        .attr("x", width/2)
        .attr("y", margin.bottom - 10)
        .style("text-anchor", "middle")
        .text(xLabel);
  
    svg.append("g")
        .classed("y axis", true)
        .call(yAxis)
      .append("text")
        .classed("label", true)
        .attr("transform", "rotate(-90)")
        .attr("x",0 - (height / 2))
        .attr("y", -margin.left +10)
        .attr("dy", ".71em")
        .style("text-anchor", "middle")
        .text(yLabel);
  
    var objects = svg.append("svg")
        .classed("objects", true)
        .attr("width", width)
        .attr("height", height);
  
    objects.append("svg:line")
        .classed("axisLine hAxisLine", true)
        .attr("x1", 0)
        .attr("y1", 0)
        .attr("x2", width)
        .attr("y2", 0)
        .attr("transform", "translate(0," + height + ")");

    objects.append("svg:line")
        .classed("axisLine vAxisLine", true)
        .attr("x1", 0)
        .attr("y1", 0)
        .attr("x2", 0)
        .attr("y2", height);
  
    objects.selectAll(".dot")
        .data(data)
      .enter().append("circle")
        .classed("dot", true)
        //.attr("r", function (d) { return 6 * Math.sqrt(d[rCat] / Math.PI); })
        .attr("r", function (d) { return 2; })
        .attr("transform", transform)
  //    .style("fill", function(d) { return color(d[colorCat]); })
        .style("fill", "steelblue")
        .on("mouseover", tip.show)
        .on("mouseout", tip.hide);
  

    svg.append("svg:text")
	.attr("class", "title")
	.attr("x", width - title.length*8)
	.attr("y", -4)
        .attr("fill", "steelblue")
        .style("font-weight","bold")
	.text(title);

    function zoom() {
      svg.select(".x.axis").call(xAxis);
      svg.select(".y.axis").call(yAxis);
  

var tickArr      = x.ticks(),
tickDistance = x(tickArr[tickArr.length - 1]) - x(tickArr[tickArr.length - 2]);

      svg.selectAll(".dot")
          .attr("transform", transform);
    }
  
    function transform(d) {
      var tickArr      = x.ticks(),
      tickDistance = x(tickArr[tickArr.length - 1]) - x(tickArr[tickArr.length - 2]);
      var pixlen = x(tickArr[1]) - x(tickArr[1]-1);
      var adj = maxBP[d[xCat]] / 2;
      var xpos =  x(d[xCat]) + (pixlen / maxBP[d[xCat]]) * (d.BP - adj);
      return "translate(" + xpos + "," + y(d[yCat]) + ")";
    }
  });
}
