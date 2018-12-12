
using System;
using System.ComponentModel;
using System.Globalization;
using System.Threading;
using System.Windows;
using System.Windows.Threading;


namespace InterfaceCrack
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly CrackSolution _resolver;
        private readonly MainViewModel _viewModel;
        private readonly BackgroundWorker _backgroundWorker;
        private NumberFormatInfo format = new NumberFormatInfo();

        private int n;
        private double a, h, nu1, nu2, mu1, mu2;
        public MainWindow()
        {
            format.NumberGroupSeparator = ",";
            format.NumberDecimalSeparator = ".";
            InitializeComponent();
            _viewModel = new MainViewModel();
            DataContext = _viewModel;
            _backgroundWorker = (BackgroundWorker)this.FindResource("BackgroundWorker");
            _resolver = new CrackSolution(PrintText, _backgroundWorker);

        }

        void PrintText(string message)
        {
            Dispatcher.BeginInvoke(DispatcherPriority.Normal, (ThreadStart)delegate 
            {
                textBox.Text += message;
                textBox.CaretIndex = textBox.Text.Length - 1;
                textBox.LineDown();
            });
        }

        private void butCalculate_Click(object sender, RoutedEventArgs e)
        {
            butCalculate.IsEnabled = false;
            n = Int32.Parse(textBoxN.Text);
            h = Double.Parse(textBoxH.Text, format);
            a = Double.Parse(textBoxA.Text, format);
            nu1 = Double.Parse(textBoxv1.Text, format);
            nu2 = Double.Parse(textBoxv2.Text, format);
            mu1 = Double.Parse(textBoxmu1.Text, format);
            mu2 = Double.Parse(textBoxmu2.Text, format);
            _backgroundWorker.RunWorkerAsync();
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            h = Double.Parse(textBoxH.Text, format);
            a = Double.Parse(textBoxA.Text, format);
            nu1 = Double.Parse(textBoxv1.Text, format);
            nu2 = Double.Parse(textBoxv2.Text, format);
            mu1 = Double.Parse(textBoxmu1.Text, format);
            mu2 = Double.Parse(textBoxmu2.Text, format);
            _viewModel.InvalidatePlot(_resolver.CalculateWithNeuralNetwork(h, a, nu1, nu2, mu1, mu2), $"h = {h}");
            _viewModel.Model.InvalidatePlot(true); 
        }

        private void backgroundWorker_DoWork(object sender, DoWorkEventArgs e)
        {
            _viewModel.InvalidatePlot(_resolver.Calculate(n, h, a, nu1, nu2, mu1, mu2), $"h = {h}");
        }

        private void backgroundWorker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            progress.Value = e.ProgressPercentage;
        }

        private void Button_Click_Clear(object sender, RoutedEventArgs e)
        {
            textBox.Text = $"Log:{Environment.NewLine}";
            _viewModel.PlotStyle = -1;
            _viewModel.Model.Series.Clear();
            _viewModel.Model.InvalidatePlot(true);
        }

        private void backgroundWorker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            _viewModel.Model.InvalidatePlot(true);
            butCalculate.IsEnabled = true;
        }
    }
}
