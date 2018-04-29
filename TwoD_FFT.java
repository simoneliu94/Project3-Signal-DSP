import java.awt.Color;

public class TwoD_FFT {
	
	final static int N = 512;
	final static Color white = new Color(255, 255, 255);
	final static Color black = new Color(0, 0, 0);
	
	public TwoD_FFT() {
		
	}
	
	public static void testSignal() {
		//Question 7a
		Picture picSignal = new Picture(N, N);

		for (int i = 220; i < 220 + 110; i++) {
			for (int j = 180; j < 180 + 140; j++) {
				picSignal.set(i, j, white);
			}
		}
	
		for (int i = 300; i < 300 + 30; i++) {
			for (int j = 205; j < 205 + 90; j++) {
				picSignal.set(i, j, black);
			}
		}
		picSignal.show();		
		
	}
	
	public static void testPulse() {
		//Question 7b
		Picture picPulse = new Picture(N, N);

		for (int i = 0; i < 30; i++) {
			for (int j = 0; j < 120; j++) {
				picPulse.set(i, j, white);
			}
		}
		for (int i = 15; i <15+15 ; i++) {
			for (int j = 15; j < 15+90; j++) {
				picPulse.set(i, j, black);
			}
		}		
		picPulse.show();
	}
	
}
