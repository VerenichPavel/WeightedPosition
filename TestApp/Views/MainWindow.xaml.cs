using Mapsui.Tiling;
using Mapsui;
using System.Windows;
using TestApp.ViewModels;

namespace TestApp.Views
{
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

            var vm = new MainViewModel();
            DataContext = vm;

            // Инициализация карты
            MapCtrl.Map = new Map();
            MapCtrl.Map.Layers.Add(OpenStreetMap.CreateTileLayer());
        }
    }
}
