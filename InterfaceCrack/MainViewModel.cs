using System;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;

namespace InterfaceCrack
{
    /// <summary>
    /// Represents the view-model for the main window.
    /// </summary>
    public class MainViewModel
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="MainViewModel" /> class.
        /// </summary>
        public MainViewModel()
        {
            // Create the plot model
            var tmp = new PlotModel { Title = "Crack representation", Subtitle = "Plot" };
            tmp.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Normalized distance y", Minimum = 0, Maximum = 1});
            tmp.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "Normalized distance x/a" });
            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            new LinearAxis().Minimum = 0;
            this.Model = tmp;
        }

        public int PlotStyle { get; set; } = -1;

        public void InvalidatePlot(double[] answers, string label)
        {
            PlotStyle++;
            //Model.Series.Clear();
            int nextHalf = 0;
            //// Create two line series (markers are hidden by default)
            //if (foolSystem)
            //{
            //    var seriesDown = new LineSeries
            //    {
            //        Title = "First half of the discrete equations",
            //        MarkerType = MarkerType.Circle
            //    };
            //    for (int i = 0; i < answers.Length / 2; i++)
            //        seriesDown.Points.Add(new DataPoint(i, answers[i]));
            //    // Add the series to the plot model
            //    Model.Series.Add(seriesDown);
            //    nextHalf = answers.Length / 2;, BrokenLineColor = OxyColors.Bisque
            //}

            var seriesUp = new LineSeries { Title = label, MarkerType = (MarkerType)PlotStyle, BrokenLineColor = OxyColor.FromUInt32((uint)new Random().Next(0,600))};
            for (int i = nextHalf; i < answers.Length; i++)
                seriesUp.Points.Add(new DataPoint(i - answers.Length / 2, 1 / (1 + Math.Exp(-answers[i]))));

            // Add the series to the plot model
            Model.Series.Add(seriesUp);
            // Axes are created automatically if they are not defined
        }

        /// <summary>
        /// Gets the plot model.
        /// </summary>
        public PlotModel Model { get; }
    }
}

