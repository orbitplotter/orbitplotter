function plotConserved(values, label) {
  if (!values || values.length === 0) {
    Plotly.newPlot("conservedPlot", []);
    return;
  }

  const trace = {
    x: [...Array(values.length).keys()],
    y: values,
    mode: "lines",
    type: "scatter",
    line: { width: 3 },
    name: label
  };

  const layout = {
    title: label,
    paper_bgcolor: "#0d0d0d",
    plot_bgcolor: "#0d0d0d",
    font: { color: "#ddd" },
    margin: { t: 30, l: 40, r: 10, b: 30 }
  };

  Plotly.newPlot("conservedPlot", [trace], layout, { displayModeBar: false });
}
