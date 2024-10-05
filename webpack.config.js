const path = require('path');

module.exports = {
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'micos.bundle.js',
    library: 'micos',
    libraryTarget: 'umd'
  },
  mode: 'development',  // Use 'development' mode for easier debugging
  devServer: {
    static: path.join(__dirname, 'src'),  // Serve files from the 'src' folder during development
    compress: true,
    port: 9000,  // Port to serve the app
    open: true,  // Automatically open in the browser
    hot: true,   // Enable hot module replacement (live reload)

  },

  module: {
    rules: [
      {
        test: /\.js$/, // Apply Babel loader to JS files        use: ['style-loader', 'css-loader'],
        exclude: /node_modules/,
        use: 'babel-loader'
      },
      
      {
        test: /\.css$/i,      // CSS rule to handle .css files
        use: ['style-loader', 'css-loader'],
      },
    ]
  }
};





