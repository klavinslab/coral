// Full JSON description of plasmid - in the future, this would be
// downloaded dynamically
function d3sequence(sequence_json, div_id) {
  // TODO: Have methods for drawing circular vs. linear, keep core logic the same
  // Option that might make it into a configurable version
  // TODO: control more of these with CSS

  var circularConfig = {
    'width': 500,
    'height': 500,
    'padding': 75,
    'featureScale': 1.2,
    'featureOpacity': 0.5,
    'backboneThickness': 0.2,
    'backboneColor': '#eee',
    'minFeatureSize': 5,
    'ticks': 80,
    'majorTickInterval': 5,
    'tickColorMinor': '#94aab0',
    'tickColorMajor': 'black',
    'tickLength': 6,
    'lineFeatureWidth': 4,
    'titleFontSize': '25px'
  };

  var linearConfig = {
    'width': 800,
    'height': 150,
    'padding': 35,
    'featureScale': 1.2,
    'featureOpacity': 0.5,
    'backboneThickness': 0.2,
    'backboneColor': '#eee',
    'minFeatureSize': 5,
    'ticks': 40,
    'majorTickInterval': 5,
    'tickColorMinor': '#94aab0',
    'tickColorMajor': 'black',
    'tickLength': 6,
    'lineFeatureWidth': 4,
    'titleFontSize': '25px'
  };

  if (!sequence_json.circular) {
    config = linearConfig;
  } else if (sequence_json.circular) {
    config = circularConfig;
  }

  ///////////////////////
  // Initialize svg space
  ///////////////////////
  var svg = d3.select('#' + div_id)
    .append('svg')
      .attr('width', config.width)
      .attr('height', config.height)
    .append('g');

  var center = {'x': config.width / 2, 'y': config.height / 2};

  if (sequence_json.circular) {

    var outerRadius = Math.min(config.width, config.height) / 2 - config.padding;
    var innerRadius = outerRadius - config.backboneThickness * outerRadius;
    var centerRadius = d3.mean([outerRadius, innerRadius]);
    var featureOffset = (config.featureScale * (outerRadius - innerRadius) - (outerRadius - innerRadius)) / 2;

    /////////////////////////////////////////
    // Placement rules for circular sequences
    /////////////////////////////////////////
    var centerTranslation = 'translate(' + center.x + ',' + center.y + ')';
    svg.attr('transform', centerTranslation);

    var baseToRadians = d3.scale.linear()
        .domain([0, sequence_json.sequence.length])
        .range([-(Math.PI / 2), 1.5 * Math.PI]);

    // Placement around the plasmid is always prior to centering in the svg
    function polarToCartesian(radius, angle) {
      return {'x': radius * Math.cos(angle), 'y': radius * Math.sin(angle)};
    }


    ///////////////////////////////////////////
    // Make an annulus for the plasmid backbone
    ///////////////////////////////////////////

    var annulus = d3.svg.arc()
        .innerRadius(innerRadius)
        .outerRadius(outerRadius)
        .startAngle(0)
        .endAngle(2 * Math.PI);

    backbone = svg.append('path')
        .attr('class', 'backbone')
        .attr('d', annulus)
        .style('fill-rule', 'evenodd')
        .style('fill', config.backboneColor)
        .style('stroke', d3.rgb(config.backboneColor).darker(0.5))
        .style('stroke-width', '1.5px');

    /////////////////////////////////
    // Add plasmid name, number of bp
    /////////////////////////////////
    svg.append('text')
        .text(sequence_json.name)
        .attr('x', 0)
        .attr('y', 0)
        .style('text-anchor', 'middle')
        .style('font-size', config.titleFontSize);

    svg.append('text')
        .text(sequence_json.sequence.length + ' bp')
        .attr('x', 0)
        .attr('y', 40)
        .style('text-anchor', 'middle');

    /////////////
    // Tick Marks
    /////////////

    // Generate ticks data
    var tickData = [];
    for (var i=0; i < config.ticks; i++) {
      tickData.push(i * Math.floor(sequence_json.sequence.length / config.ticks));
    }

    // FIXME: why use paths to make Ticks when lines can be used?
    var line = d3.svg.line()
        .x(function(d) { return d.x; } )
        .y(function(d) { return d.y; } );

    tickGroup = svg.append('g')
        .attr('id', 'tick-group');

    ticks = tickGroup.selectAll('.tick')
        .data(tickData)
      .enter().append('path')
        .attr('d', function(d, i) {
          tickAngle = baseToRadians(d);
          inner = polarToCartesian(outerRadius, tickAngle);
          outer = polarToCartesian(outerRadius + config.tickLength, tickAngle);
          return line([inner, outer]);
        })
        .attr('stroke', function(d, i) {
          if (i === 0) {
            tickColor = 'red';
          } else if (i % config.majorTickInterval === 0) {
            tickColor = config.tickColorMajor;
          } else {
            tickColor = config.tickColorMinor;
          }
          return tickColor;
        })
        .attr('stroke-width', 2)
        .attr('fill', 'none');

    tickLabels = tickGroup.selectAll('.tick-label')
        .data(tickData)
      .enter().append('text')
        .text(function(d) {
            if (d % config.majorTickInterval === 0) {
              return d;
            }
        })
        .attr('x', function(d) {
          tickAngle = baseToRadians(d);
          return polarToCartesian(outerRadius + 22, tickAngle).x;
        })
        .attr('y', function(d) {
          tickAngle = baseToRadians(d);
          return polarToCartesian(outerRadius + 22, tickAngle).y;
        })
        .attr('dy', '0.4em')
        .attr('fill', config.tickColorMinor)
        .style('text-anchor', 'middle')
        .style('font-size', '12px');

    ////////////////////////////////
    // Add features and their labels
    ////////////////////////////////

    // Separate out arc features and line features
    var blockData = [];
    var lineData = [];
    for (var i = 0; i < sequence_json.features.length; i++) {
      feature = sequence_json.features[i];
      arcAngle = 2 * Math.PI / sequence_json.sequence.length * (feature.stop - feature.start);
      arcLength = centerRadius * arcAngle;
      if (arcLength > config.minFeatureSize) {
        blockData.push(feature);
      } else {
        lineData.push(feature);
      }
    }

    var featuresGroup = svg.append('g')
        .attr('id', 'features');
    var blockGroup = featuresGroup.append('g')
        .attr('id', 'block-features');

    var blockFeatures = blockGroup.selectAll('.feature')
        .data(blockData)
      .enter()
        .append('g')
        .attr('class', 'feature')
        .append('path')
        .attr('d', function(d, i) {
          the_arc = d3.svg.arc()
            // d3 svg arc origin is at the top
            .startAngle(baseToRadians(d.start) + Math.PI / 2)
            .endAngle(baseToRadians(d.stop) + Math.PI / 2)
            .innerRadius(innerRadius - featureOffset)
            .outerRadius(outerRadius + featureOffset)
            (d, i);
          return the_arc;
        })
        .attr('id', function(d, i) { return 'block' + i; })
        .style('fill-rule', 'evenodd')
        .style('fill', function(d) {
          if (d.hasOwnProperty('color')) {
            return d.color;
          } else {
            return 'rgb(192, 192, 256)';
          }
        })
        .style('fill-opacity', config.featureOpacity)
        .style('stroke', function(d) {
          if (d.hasOwnProperty('color')) {
            color = d.color;
          } else {
            color = 'rgb(192, 192, 256)';
          }
          return d3.rgb(color).darker(0.5);
        })
        .style('stroke-width', '1.5px');


    // Gaussian blur filter for making text shadow
    var filter = svg.append('defs')
      .append('filter')
        .attr('id', 'text-shadow')
        .attr('height', '130%')
      .append('feGaussianBlur')
        .attr('in', 'SourceAlpha')
        .attr('stdDeviation', 2)
        .attr('result', 'blur');

    // Add invisible arcs for aligning text
    function describeArc(radius, start, end) {
      start_xy = polarToCartesian(radius, baseToRadians(start));
      end_xy = polarToCartesian(radius, baseToRadians(end));

      var d = [
        'M', start_xy.x, start_xy.y,
        'A', radius, radius, 0, 0, 1, end_xy.x, end_xy.y
      ].join(' ');

      return d;
    };

    var centerArcs = blockGroup.selectAll('.feature, .center-arcs')
      .append('path')
        .attr('d', function(d) {
          return describeArc(0.97 * centerRadius, d.start, d.stop)
        })
        .attr('class', 'center-arc')
        .attr('id', function(d, i) { return div_id + '-center-arc' + i; })
        .style('stroke', 'none')
        .style('fill', 'none');

    var shadows = blockGroup.selectAll('.feature')
      .append('text')
        .attr('class', 'block-label-shadow')
      .append('textPath')
        .attr('class', 'shadow')
        .attr('startOffset', '50%')
        .attr('xlink:href', function(d, i) {
          return '#' + div_id + '-center-arc' + i;
        })
        .attr('id', function(d, i) {
          return 'shadow' + i;
        })
        .style('text-anchor', 'middle')
        .style('font-family', 'sans-serif')
        .style('font-weight', 'bold')
        .style('font-size', '11px')
        .style('filter', 'url(#text-shadow)')
        .text(function(d) { return d.name; });

    var blockText = blockGroup.selectAll('.feature')
      .append('text')
        .attr('class', 'block-label')
      .append('textPath')
        .attr('class', 'block-text')
        .attr('fill', 'white')
        .attr('startOffset', '50%')
        .attr('xlink:href', function(d, i) {
          return '#' + div_id + '-center-arc' + i;
        })
        .attr('id', function(d, i) { return 'label' + i; })
        .style('text-anchor', 'middle')
        .style('font-family', 'sans-serif')
        .style('font-weight', 'bold')
        .style('font-size', '11px')
        .text(function(d) { return d.name; });

    // If text is too big, don't offset and set anchor to left
    // FIXME: this should really be done with more d3-style methods
    feature_selection = document.getElementsByClassName('feature');
    for (var i = 0; i < feature_selection.length; i++) {
      blocktemp = feature_selection[i];
      text = blocktemp.getElementsByClassName('block-text')[0];
      shadow = blocktemp.getElementsByClassName('shadow')[0];

      centerArc = blocktemp.getElementsByClassName('center-arc')[0];
      blockLength = centerArc.getTotalLength();

      nameLength = text.innerHTML.length;
      while (1.2 * text.getComputedTextLength() > blockLength) {
         nameLength -= 1;
         // Trim until it isn't
         text.innerHTML = text.innerHTML.slice(0, nameLength);
         text.setAttribute('startOffset', '10%');
         text.style['text-anchor'] = 'start';

         shadow.innerHTML = shadow.innerHTML.slice(0, nameLength);
         shadow.setAttribute('startOffset', '10%');
         shadow.style['text-anchor'] = 'start';
       }
    }

    // Add line features
    var lineGroup = featuresGroup.append('g')
        .attr('class', 'lineFeatureGroup');

    var lineFeatures = lineGroup.selectAll('.feature')
        .data(lineData)
      .enter().append('g')
      .append('path')
        .attr('class', 'feature')
        .attr('class', 'line-feature')
        .attr('d', function(d, i) {
          featureAngle = baseToRadians((d.start + d.stop) / 2);
          inner = polarToCartesian(innerRadius - featureOffset - config.tickLength, featureAngle);
          outer = polarToCartesian(outerRadius + featureOffset + config.tickLength, featureAngle);

          return line([inner, outer]);
        })
        .attr('stroke', 'rgb(96, 96, 96)')
        .attr('stroke-width', config.lineFeatureWidth)
        .attr('fill', 'none')
        .on('mouseover', function(d) {
          angle = baseToRadians((d.start + d.stop) / 2);
          position = polarToCartesian((outerRadius + featureOffset + config.tickLength + 3), angle);
          lineGroup.append('text')
              .attr('id', 'lineFeatureHighlight')
              .attr('x', position.x)
              .attr('y', position.y)
              .attr('text-anchor', 'middle')
              .attr('font-family', 'sans-serif')
              .attr('font-size', '14px')
              .style('font-weight', 'bold')
              .attr('fill', 'red')
              .attr('dy', '-1em')
              .attr('transform', function(d) {
                degrees = (angle + Math.PI) / (2 * Math.PI) * 360;
                return 'rotate(' + (degrees - 90) + ' ' + (position.x) + ',' + (position.y) + ')';
              })
              .text(d.name);
        })
        .on('mouseout', function(d) {
          d3.select('#lineFeatureHighlight').remove();
        });
    } else {

    ///////////////////////////////////////
    // Placement rules for linear sequences
    ///////////////////////////////////////
    var backboneHeight = 25;  // TODO: replace with unified thickness param
    var xMin = config.padding;
    var xMax = config.width - config.padding;

    var backboneDim = {
      'left': xMin,
      'right': xMax,
      'top': center.y - backboneHeight / 2,
      'bottom': center.y + backboneHeight / 2
    };

    // Offset of 1 pixel makes the display look nicer
    var xScale = d3.scale.linear()
        .domain([0, sequence_json.sequence.length])
        .range([xMin + 1, xMax + 1]);

    ////////////////////
    // Draw the backbone
    ////////////////////
    // TODO: similarities with circular - same fill, stroke, etc. Could simplify
    //       with a backboneWidth attribute that's in pixels
    svg.append('rect')
        .attr('class', 'backbone')
        .attr('x', backboneDim.left)
        .attr('y', backboneDim.top)
        .attr('width', backboneDim.right - backboneDim.left)
        .attr('height', backboneDim.bottom - backboneDim.top)
        .style('fill', config.backboneColor)
        .style('stroke', d3.rgb(config.backboneColor).darker(0.5))
        .style('stroke-width', '1.5px');

    /////////////////////////////////
    // Add plasmid name, number of bp
    /////////////////////////////////
    // TODO: similarities with circular - identical, just different placement
    // FIXME: make the placement of text have no magic numbers
    svg.append('text')
        .text(sequence_json.name)
        .attr('x', center.x)
        .attr('y', backboneDim.bottom +  50)
        .style('text-anchor', 'middle')
        .style('font-size', config.titleFontSize);

    svg.append('text')
        .text(sequence_json.sequence.length + ' bp')
        .attr('x', center.x)
        .attr('y', backboneDim.bottom + 80)
        .style('text-anchor', 'middle');

    /////////////
    // Tick Marks
    /////////////
    // TODO: similarities with circular - almost everything
    // Generate ticks data
    var tickData = [];
    for (var i=0; i < config.ticks; i++) {
      tickData.push(i * Math.floor(sequence_json.sequence.length / config.ticks));
    }

    var tickGroup = svg.append('g')
        .attr('id', 'tick-group');

    var ticks = tickGroup.selectAll('.tick')
        .data(tickData)
      .enter().append('line')
        .attr('x1', function(d) { return xScale(d); })
        .attr('x2', function(d) { return xScale(d); })
        .attr('y1', backboneDim.bottom)
        .attr('y2', backboneDim.bottom + config.tickLength)
        .attr('stroke', function(d, i) {
          if (i % config.majorTickInterval === 0) {
            tickColor = config.tickColorMajor;
          } else {
            tickColor = config.tickColorMinor;
          }
          return tickColor;
        })
        .attr('stroke-width', 2)
        .attr('fill', 'none');

    var tickLabels = tickGroup.selectAll('.tick-label')
        .data(tickData)
      .enter().append('text')
        .text(function(d) {
            if (d % config.majorTickInterval === 0) {
              return d;
            }
        })
        .attr('x', function(d) { return xScale(d); })
        .attr('y', function(d) { return backboneDim.bottom + config.tickLength + 8; })
        .attr('dy', '0.4em')
        .attr('fill', config.tickColorMinor)
        .style('text-anchor', 'middle')
        .style('font-size', '12px');

    ////////////////////////////////
    // Add features and their labels
    ////////////////////////////////

    // Separate out arc features and line features
    var blockData = [];
    var lineData = [];
    // FIXME: this process is broken and should be determined by fit, not n
    // arbitrary scaling factor
    for (var i = 0; i < sequence_json.features.length; i++) {
      feature = sequence_json.features[i];
      if (xScale(feature.stop - feature.start) > 10 * config.minFeatureSize) {
        blockData.push(feature);
      } else {
        lineData.push(feature);
      }
    }

    var featuresGroup = svg.append('g')
        .attr('id', 'features');
    var blockGroup = featuresGroup.append('g')
        .attr('id', 'block-features');

    var blockFeatures = blockGroup.selectAll('.feature')
        .data(blockData)
      .enter().append('g')
        .attr('class', 'feature')
      .append('rect')
        .attr('class', 'feature-block')
        .attr('x', function(d) { return xScale(d.start); })
        .attr('y', function(d) { return backboneDim.top - (config.featureScale * backboneHeight - backboneHeight) / 2; })
        .attr('width', function(d) { return xScale(d.stop) - xScale(d.start); })
        .attr('height', function(d) { return config.featureScale * backboneHeight; })
        .attr('id', function(d, i) { return 'block' + i; })
        .style('fill', function(d) {
          if (d.hasOwnProperty('color')) {
            return d.color;
          } else {
            return 'rgb(192, 192, 256)';
          }
        })
        .style('fill-opacity', config.featureOpacity)
        .style('stroke', function(d) {
          if (d.hasOwnProperty('color')) {
            color = d.color;
          } else {
            color = 'rgb(192, 192, 256)';
          }
          return d3.rgb(color).darker(0.5);
        })
        .style('stroke-width', '1.5px');

    // Gaussian blur filter for making text shadow
    var filter = svg.append('defs')
      .append('filter')
        .attr('id', 'text-shadow')
        .attr('height', '130%')
      .append('feGaussianBlur')
        .attr('in', 'SourceAlpha')
        .attr('stdDeviation', 2)
        .attr('result', 'blur');

    var shadows = blockGroup.selectAll('.feature')
      .append('text')
        .attr('class', 'shadow')
        .attr('id', function(d, i) {
          return 'shadow' + i;
        })
        .attr('x', function(d) { return xScale((d.start + d.stop) / 2); })
        .attr('y', function(d) { return (backboneDim.top + backboneDim.bottom) / 2 + 3; })
        .style('text-anchor', 'middle')
        .style('font-family', 'sans-serif')
        .style('font-weight', 'bold')
        .style('font-size', '10px')
        .style('filter', 'url(#text-shadow)')
        .text(function(d) { return d.name; });

    var blockText = blockGroup.selectAll('.feature')
      .append('text')
        .attr('class', 'block-text')
        .attr('fill', 'white')
        .attr('id', function(d, i) { return 'label' + i; })
        .attr('x', function(d) { return xScale((d.start + d.stop) / 2); })
        .attr('y', function(d) { return (backboneDim.top + backboneDim.bottom) / 2; })
        .style('text-anchor', 'middle')
        .style('font-family', 'sans-serif')
        .style('font-weight', 'bold')
        .style('font-size', '10px')
        .text(function(d) { return d.name; });

    // If text is too big, don't offset and set anchor to left
    // FIXME: this should really be done with more d3-style methods
    feature_selection = document.getElementsByClassName('feature');
    for (var i = 0; i < feature_selection.length; i++) {
      blocktemp = feature_selection[i];
      text = blocktemp.getElementsByClassName('block-text')[0];
      shadow = blocktemp.getElementsByClassName('shadow')[0];

      block = blocktemp.getElementsByClassName('feature-block');
      blockLength = block[0].width.baseVal.value

      nameLength = text.innerHTML.length;
      while (1.2 * text.getComputedTextLength() > blockLength) {
         nameLength -= 1;
         // Trim until it isn't
         text.innerHTML = text.innerHTML.slice(0, nameLength);

         shadow.innerHTML = shadow.innerHTML.slice(0, nameLength);
       }
    }

    // Add line features
    var lineGroup = featuresGroup.append('g')
        .attr('class', 'lineFeatureGroup');

    var lineFeatures = lineGroup.selectAll('.feature')
        .data(lineData)
      .enter().append('g')
      .append('line')
        .attr('class', 'feature')
        .attr('class', 'line-feature')
        .attr('x1', function(d) { return xScale(d.start); })
        .attr('x2', function(d) { return xScale(d.start); })
        .attr('y1', function(d) { return backboneDim.bottom + (config.featureScale * backboneHeight - backboneHeight) / 2 + config.tickLength; })
        .attr('y2', function(d) { return backboneDim.top - (config.featureScale * backboneHeight - backboneHeight) / 2 - config.tickLength; })
        .attr('stroke', 'rgb(96, 96, 96)')
        .attr('stroke-width', config.lineFeatureWidth)
        .attr('fill', 'none')
        .on('mouseover', function(d) {
          lineGroup.append('text')
              .attr('id', 'lineFeatureHighlight')
              .attr('x', xScale(d.start))
              .attr('y', backboneDim.top - (config.featureScale * backboneHeight - backboneHeight) / 2 - 3)
              .attr('text-anchor', 'middle')
              .attr('font-family', 'sans-serif')
              .attr('font-size', '14px')
              .style('font-weight', 'bold')
              .attr('fill', 'red')
              .attr('dy', '-1em')
              .text(d.name);
        })
        .on('mouseout', function(d) {
          d3.select('#lineFeatureHighlight').remove();
        });

  }
}
