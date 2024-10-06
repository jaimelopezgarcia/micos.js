const path = require('path');

module.exports = {
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'micos.bundle.js',
    library: 'micos',
    libraryTarget: 'umd',
    publicPath: './', // Ensure relative paths
  },
  mode: 'development',
  devServer: {
    static: path.join(__dirname, '/'),
    compress: true,
    port: 9000,
    open: true,
    hot: true,
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        exclude: /node_modules/,
        use: 'babel-loader',
      },
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
    ],
  },
};
