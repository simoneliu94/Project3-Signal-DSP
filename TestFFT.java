
public class TestFFT {
	public static void main(String[] args){
		double[][] data = {{26160.0},
				{19011.0},
				{18757.0},
				{18405.0},
				{17888.0},
				{14720.0},
				{14285.0},
				{17018.0},
				{18014.0},
				{17119.0},
				{16400.0},
				{17497.0},
				{17846.0},
				{15700.0},
				{17636.0},
				{17181.0}};
				
		Matrix mat = new Matrix(16, 1);
		mat.matrix = data;
		
		System.out.println("Data");
		System.out.println(mat);
		
		System.out.println("FFT of Data");
		mat.fft(1);				
		for(int i = 0; i < mat.numRows; i++){
			System.out.println(mat.complexMatrix[i]);
		}
				
	}
}
