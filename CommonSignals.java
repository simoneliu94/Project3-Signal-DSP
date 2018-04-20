import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

public class CommonSignals {
	final static int SAMPLES = 512;
	
	
	public static Matrix fs(int s){
		Matrix f = new Matrix(SAMPLES, 1);
		double f_total = 0;
		double t = 0;
		
		f.matrix[0][0] = 0;
		for(int j=0; j<SAMPLES-1; j++){	
			f_total = 0;
			t = t + .001953125; // 1/SAMPLES
			for (int k=1; k<=s; k++){
				f_total = f_total + (Math.sin(2*Math.PI*(2*k-1)*t)/(2*k-1));
			}		

			f.matrix[j+1][0] = f_total;
			
		}
		//System.out.println(f);
		return f;
	}
	
	
	public static Matrix gs(int s){
		Matrix g = new Matrix(SAMPLES, 1);
		double g_total = 0;
		double t = 0;		
		
		g.matrix[0][0] = 0;
		for(int j=0; j<SAMPLES-1; j++){	
			g_total = 0;
			t = t + .001953125; // 1/SAMPLES
			for (int k=1; k<=s; k++){
				g_total = g_total + (Math.sin(2*Math.PI*(2*k)*t)/(2*k));
			}		

			g.matrix[j+1][0] = g_total;
			
		}
		return g;
		
	}
	
	public static Matrix generate_FL(int s) {
		Matrix f = new Matrix (SAMPLES, 1);
		double t = 0;
		double data = 0;
		
		f.matrix[0][0] = 0;
		for(int i=0; i<SAMPLES-1; i++){	
			t = t + .001953125;
			if (t>0 && t<0.5){
				data = Math.PI/4;
			}
			
			else if (t>0.5 && t<1){
				data = -Math.PI/4;
			}
			
			else if(t == 0.0 || t == 0.5){
				data = 0; 
			}
			f.matrix[i+1][0] = data;
		}
		return f;
	}
	
	public static Matrix generate_GL(int s) {
		Matrix g = new Matrix (SAMPLES, 1);
		double t = 0;
		double data = 0;
		
		g.matrix[0][0] = 0;
		for(int i=0; i<SAMPLES-1; i++){	
			t = t + .001953125;
			if (t>0 && t<0.5){
				data = -(Math.PI)*t + (Math.PI/4);
			}
			
			else if (t>0.5 && t<1){
				data = -(Math.PI)*(t-0.5) + (Math.PI/4);
			}
			
			else if(t == 0.0 || t == 0.5){
				data = 0; 
			}
			g.matrix[i+1][0] = data;
		}
		return g;
	}
	
	public static void printComplex(Complex[] cList) {
		for(int i = 0; i < cList.length; i++){
			System.out.println(cList[i]);
		}
	}
	
	public static void printComplex_real(Complex[] cList) {
		for(int i = 0; i < cList.length; i++){
			System.out.println(cList[i].re);
		}
	}
	
	public static Complex[] lowpass(Complex[] f_fft) {
		Matrix lowpass = new Matrix(512,1);
		for(int i = 0; i < lowpass.numRows; i++){
			lowpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i = 0; i < 14; i++){
			lowpass.complexMatrix[i] = new Complex(1,0);
		}
		for(int i = lowpass.numRows - 14; i < lowpass.numRows; i++){ //7 or 14
			lowpass.complexMatrix[i] = new Complex(1,0);
		}		
		
		for(int i = 0; i < lowpass.numRows; i ++){
			lowpass.complexMatrix[i] = lowpass.complexMatrix[i].times(f_fft[i]);
		}		
		
		return lowpass.fft(-1);
	}
	
	public static Complex[] highpass(Complex[] f_fft) {
		Matrix highpass = new Matrix(512,1);
		for(int i = 0; i < highpass.numRows; i++){
			highpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i = 0; i < 14; i++){
			highpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i = highpass.numRows- 14; i < highpass.numRows; i++){
			highpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i = 14; i<100; i++) {
			highpass.complexMatrix[i] = new Complex(1,0);
		}
		for(int i = 100; i < highpass.numRows; i++){
			highpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i=highpass.numRows-100; i<highpass.numRows-14; i++) {
			highpass.complexMatrix[i] = new Complex(1,0);
		}
		
		for(int i = 0; i < highpass.numRows; i ++){
			highpass.complexMatrix[i] = highpass.complexMatrix[i].times(f_fft[i]);
		}			
		return highpass.fft(-1);
	}
	
	public static Complex[] bandpass(Complex[] f_fft) {
		Matrix bandpass = new Matrix(512,1);
		for(int i = 0; i < bandpass.numRows; i++){
			bandpass.complexMatrix[i] = new Complex(0,0);
		}
		for(int i = 0; i < bandpass.numRows; i++){//11,13,15,17
			if (i==9 || i==11 ||i==13 ||i==15 ||i==497 ||i==499 ||i==501 ||i==503) {
				bandpass.complexMatrix[i] = new Complex(1,0);
			}			
		}		
		for(int i = 0; i < bandpass.numRows; i ++){
			bandpass.complexMatrix[i] = bandpass.complexMatrix[i].times(f_fft[i]);
		}			
		return bandpass.fft(-1);
	}
	
	public static Complex[] notchpass(Complex[] f_fft) {
		Matrix notchpass = new Matrix(512,1);
		for(int i = 0; i < notchpass.numRows; i++){
			notchpass.complexMatrix[i] = new Complex(1,0);
		}
		for(int i = 0; i < notchpass.numRows; i++){//11,13,15,17
			if (i==9 || i==11 ||i==13 ||i==15 ||i==497 ||i==499 ||i==501 ||i==503) {
				notchpass.complexMatrix[i] = new Complex(0,0);
			}			
		}		
		for(int i = 0; i < notchpass.numRows; i ++){
			notchpass.complexMatrix[i] = notchpass.complexMatrix[i].times(f_fft[i]);
		}			
		return notchpass.fft(-1);
	}

	
	public static void main(String[] args) throws IOException{
		//Question1
		//fs(3);
		//fs(10);
		fs(50);
		//gs(3);
		//gs(10);
		//gs(50);	
		
		//fs(7);
		//System.out.println(fs(7));
		
		//System.out.println(fs(50).subtract(fs(8)).add(fs(4)));
		
		
		//FFT of f50
		//Matrix f50 = new Matrix (fs(50));
		//Complex[] f50_fft = f50.fft(1);	
				
		/*System.out.println("FFT of f50");
		printComplex(f50_fft);
		
		//PSD of f50
		System.out.println("-----------------------------------------");
		System.out.println("PSD of f50");
		Matrix f50_psd = f50.psd();
		System.out.println(f50_psd);
		
		//--------------------------------------------------------------------------------------
		
		//FFT of g50
		Matrix g50 = new Matrix (gs(50));
		Complex[] g50_fft = g50.fft(1);	
				
		System.out.println("FFT of g50");
		printComplex(g50_fft);
		
		//PSD of f50
		System.out.println("-----------------------------------------");
		System.out.println("PSD of g50");
		Matrix g50_psd = g50.psd();
		System.out.println(g50_psd);*/
		
		/*Matrix fL = new Matrix (generate_FL(50));
		System.out.println("Data of fL(50)");
		System.out.println(fL);
		
		Matrix gL = new Matrix (generate_GL(50));
		System.out.println("Data of gL(50)");
		System.out.println(gL);		
		
		fL.fft(1);
		gL.fft(1);
		
		System.out.println("-----------------------------------------");
		System.out.println("PSD of fL50");
		Matrix fL50_psd = fL.psd();
		System.out.println(fL50_psd);
		
		System.out.println("-----------------------------------------");
		System.out.println("PSD of gL50");
		Matrix gL50_psd = gL.psd();
		System.out.println(gL50_psd);*/
		
		
		
		//Question2
		/*System.out.println("Sum of 2 sin functions");
		Matrix sinSum = new Matrix(1, "Resources/sinSum.txt");		
		Matrix sinSum_matrix = new Matrix (sinSum);
		Complex[] sinSum_fft = sinSum_matrix.fft(1);
		printComplex(sinSum_fft);
		System.out.println("-----------------------------------------");
		System.out.println("PSD of sinSum");
		Matrix sinSum_psd = sinSum_matrix.psd();
		System.out.println(sinSum_psd);*/
		
		
		/*System.out.println("Product of 2 sin functions");
		Matrix sinProd = new Matrix(1, "Resources/sinProduct.txt");		
		Matrix sinProd_matrix = new Matrix (sinProd);
		Complex[] sinProd_fft = sinProd_matrix.fft(1);
		printComplex(sinProd_fft);
		System.out.println("-----------------------------------------");
		System.out.println("PSD of sinProduct");
		Matrix sinProd_psd = sinProd_matrix.psd();
		System.out.println(sinProd_psd);*/
		
		//Question4
		//printComplex_real(lowpass(f50_fft));
		//printComplex_real(highpass(f50_fft));
		//printComplex_real(bandpass(f50_fft));
		//printComplex_real(notchpass(f50_fft));
		

	}
}
