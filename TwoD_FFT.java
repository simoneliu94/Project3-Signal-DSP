import java.awt.Color;

public class TwoD_FFT {
	public Complex[][] TwoDFFT() {

		Complex[][] pixelsSignal = new Complex[512][512];
		Complex[][] pixelsPulse = new Complex[512][512];

		int N = 512;
		Color white = new Color(255, 255, 255);

		Picture picSignal = new Picture(N, N);
		Picture picOutput = new Picture(N, N);

		for (int i = 230; i < 230 + 150; i++) {

			for (int j = 260; j < 260 + 75; j++) {

				picSignal.set(i, j, white);

			}// end j for

		}// end i for

		picSignal.show();

		Picture picPulse = new Picture(N, N);

		for (int i = 0; i < 50; i++) {

			for (int j = 0; j < 25; j++) {

				picPulse.set(i, j, white);

			}// end j for

		}// end i for

		picPulse.show();

		// get the signal colors

		Color[][] colorArray = new Color[512][512];

		colorArray = picSignal.getColorArray();

		for (int m = 0; m < colorArray.length; m++) {

			for (int l = 0; l < colorArray[0].length; l++) {

				pixelsSignal[m][l] = new Complex(colorArray[m][l].getRed(), 0);

			}// end l

		}// end m

		// get the pusle colors

		colorArray = picPulse.getColorArray();

		for (int m = 0; m < colorArray.length; m++) {

			for (int l = 0; l < colorArray[0].length; l++) {

				pixelsPulse[m][l] = new Complex(colorArray[m][l].getRed(), 0);

			}// end l

		}// end m

		/* Compute m FFTs of length n */

		Complex[][] pixelsSignalFFT = new Complex[512][512];
		Complex[][] pixelsPulseFFT = new Complex[512][512];
		Complex[][] pixelsConjugate = new Complex[512][512];
		Complex[][] productTwoD = new Complex[512][512];
		Complex[][] productInv = new Complex[512][512];

		Complex[][] pixelsPulseTemp = new Complex[512][512];

		for (int i = 0; i < pixelsPulseTemp.length; i++) {

			for (int j = 0; j < pixelsPulseTemp[0].length; j++) {

				pixelsPulseTemp[i][j] = pixelsPulse[i][j];

			}// end j for

		}// end i for
			// //////////////////////////////////////////////////////Signal////////////////////////////

		pixelsSignalFFT = algorithm(pixelsSignal, 1);
		pixelsPulseFFT = algorithm(pixelsPulse, 1);
		pixelsConjugate = conjugateTwoD(pixelsPulseFFT);

		System.out.println("Pulse");
		// printtwod(pixelsPulse);
		System.out.println("Signal");
		// printtwod(pixelsSignal);

		//
		// The conjugate of the pulseFFT

		// Multiply the signal with the pulse conjugate
		productTwoD = multiplyFFTTwoD(pixelsSignalFFT, pixelsConjugate);

		// ////////////////////////////////////////////////////////////////////////////////////
		// Take the inverse of the product
		productInv = algorithm(productTwoD, -1);

		// Find the max of the inverse

		double max = -999999999;
		double min = 999999999;

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				if (productInv[i][j].x > max) {

					max = productInv[i][j].x;

				}// end if

			}// end j for

		}// end i for

		// ////////////////////////////////ok/////////////////////

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				if (productInv[i][j].x < min) {

					min = productInv[i][j].x;

				}// end if

			}// end j
		}// end i

		printtwod(productInv);

		System.out.println("max = " + max);
		System.out.println("min = " + min);

		// double mx = Math.log(max - min + 2);
		double mx = (255 / max);
		double out;

		// ////////////////////////////////////////////////////

		for (int i = 0; i < productInv.length; i++) {

			for (int j = 0; j < productInv[0].length; j++) {

				// out = Math.log(productInv[i][j].x - min + 1) * (255 / mx);
				out = productInv[i][j].x * mx;

				productInv[i][j].setX((int) out);

			}// end j for

		}// end i for

		printtwod(productInv);

		Picture resultingPicture = new Picture(512, 512);

		for (int i = 0; i < resultingPicture.height(); i++) {

			for (int j = 0; j < resultingPicture.width(); j++) {

				Color color = new Color((int) productInv[i][j].getX(),
						(int) productInv[i][j].getX(),
						(int) productInv[i][j].getX());

				resultingPicture.set(i, j, color);

			}// end j for

		}// end i for

		resultingPicture.show();

		Picture redPicture = new Picture(512, 512);

		for (int i = 0; i < redPicture.height(); i++) {

			for (int j = 0; j < redPicture.width(); j++) {

				if (productInv[i][j].getX() >= 229.5
						&& productInv[i][j].getX() < 255) {

					Color red = new Color(255, 0, 0);
					redPicture.set(i, j, red);

				} else {

					Color color = new Color((int) productInv[i][j].getX(),
							(int) productInv[i][j].getX(),
							(int) productInv[i][j].getX());

					redPicture.set(i, j, color);

				}// end else

			}// end j for
		}// end i for

		redPicture.show();

		// printtwod(product);
		//
		// for (i = 0; i < product.length; i++) {
		//
		// for (int j = 0; j < product[0].length; j++) {
		//
		// if (product[i][j].x == 255) {
		//
		// picOutput.set(i, j, Color.WHITE);
		//
		// } else {
		//
		// picOutput.set(i, j, Color.BLACK);
		//
		// }

		// }// end j for
		//
		// }// end i for

		// System.out.println("Showing output");
		// picOutput.show();

		return null;

	}// end TwoDFFT
	
	public Complex[][] conjugateTwoD(Complex[][] a) {

		Complex[][] z = new Complex[a.length][a[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a[0].length; j++) {

				z[i][j] = a[i][j];

			}// end j for

		}// end i for

		for (int i = 0; i < z.length; i++) {

			for (int j = 0; j < z[0].length; j++) {

				z[i][j].negateY();

			}// end j for

		}// end i for

		return z;

	}// end conjugate
	
	public Complex[][] algorithm(Complex[][] pixelsArray, int b) {

		// //////////////////////////////////////////////////////Signal///////////////////////////

		int d = b;

		Complex[][] Y;
		Complex[][] Z = new Complex[pixelsArray.length][pixelsArray[0].length];

		for (int i = 0; i < pixelsArray.length; i++) {

			for (int j = 0; j < pixelsArray[0].length; j++) {

				Z[i][j] = pixelsArray[i][j];

			}// end j for

		}// end i for
		int m = pixelsArray.length, n = pixelsArray[0].length;

		/* Initialize */

		Complex[] row1FFT = new Complex[pixelsArray.length];
		Complex[] col1FFT = new Complex[pixelsArray.length];
		Complex[] row = new Complex[pixelsArray.length];
		Complex[] col = new Complex[pixelsArray.length];
		Y = new Complex[pixelsArray.length][pixelsArray[0].length];

		int i, j, k;
		/* Compute m FFTs of length n */

		/* Compute m FFTs of length n */

		for (k = 0; k < m; k++) {
			for (i = 0; i < n; i++)
				row[i] = Z[k][i];
			row1FFT = FFT(row, n, d);
			for (i = 0; i < n; i++)
				Y[i][k] = row1FFT[i];
		}

		/* Compute n FFTs of length m */

		for (i = 0; i < n; i++) {
			for (k = 0; k < m; k++)
				col[k] = Y[i][k];
			col1FFT = FFT(col, m, d);
			for (k = 0; k < m; k++)
				Y[i][k] = col1FFT[k];
		}

		/* Finalize */

		for (k = 0; k < m; k++)
			for (i = 0; i < n; i++)
				Z[k][i] = Y[i][k];

		return Z;

	}
	
	public Complex[][] multiplyFFTTwoD(Complex[][] a, Complex[][] b) {

		Complex[][] output = new Complex[a.length][a[0].length];

		Complex[][] first = new Complex[a.length][a[0].length];
		Complex[][] second = new Complex[b.length][b[0].length];

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < a.length; j++) {

				first[i][j] = a[i][j];

			}// end j for
		}// end i for

		for (int i = 0; i < b.length; i++) {

			for (int j = 0; j < b.length; j++) {

				second[i][j] = b[i][j];

			}// end j for
		}// end i fo
		if (a.length != b.length) {

			System.err.println("The two Complex arrays do not match");

		} else {

			for (int i = 0; i < a.length; i++) {

				for (int j = 0; j < a[0].length; j++) {

					second[i][j].negateY();
					output[i][j] = Complex.cmul(first[i][j], second[i][j]);

				}// end j for

			}// end i for

		}// end else

		return output;

	}// end multiplyFFTTwoD

}
