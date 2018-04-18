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
	


	
	public static void main(String[] args) throws IOException{
		//Question1
		//fs(3);
		//fs(10);
		//fs(50);
		//gs(3);
		//gs(10);
		//gs(50);		
		
		//FFT of f50
		/*Matrix f50 = new Matrix (fs(50));
		Complex[] f50_fft = f50.fft(1);	
				
		System.out.println("FFT of f50");
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
		
		
		System.out.println("Product of 2 sin functions");
		Matrix sinProd = new Matrix(1, "Resources/sinProduct.txt");		
		Matrix sinProd_matrix = new Matrix (sinProd);
		Complex[] sinProd_fft = sinProd_matrix.fft(1);
		printComplex(sinProd_fft);
		System.out.println("-----------------------------------------");
		System.out.println("PSD of sinProduct");
		Matrix sinProd_psd = sinProd_matrix.psd();
		System.out.println(sinProd_psd);
		
		

	}
}
