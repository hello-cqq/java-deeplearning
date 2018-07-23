package deeplearning.math;

import java.util.Arrays;

/**
 * 深度学习中矩阵的相关计算
 * @author cqq	1324095615@qq.com
 * @version 1.0
 * 
 */
public class Matrix {
	/**
	 * 获取n维单位矩阵
	 * @param n 指定单位矩阵的维数
	 * @return 返回n维单位矩阵
	 */
	public  static double[][] getIdentityMatrix( int n)
	{
		double[][] Identity_Matrix;
		if(n>0){
			Identity_Matrix=new double[n][n];
			for(int i=0;i<n;i++)Identity_Matrix[i][i]=1;
			return Identity_Matrix;
		}
		else {
			System.out.println("输入不合法，请确保输入参数n为大于0的int型整数");
			return null;
		}
	}
	
	/**
	 * 对单个矩阵进行简单的数学计算：与标量相乘，与标量相加
	 * @param m 输入矩阵or数组
	 * @param w 与矩阵相乘的系数
	 * @param c 与矩阵相加的系数
	 * @return 返回运算后的矩阵
	 */
	public static double[][] matrixCompute(double[][] m,double w,double c)
	{
		double[][] m_result=null;//定义返回矩阵
		try {
			m_result = new double[m.length][m[0].length];
			for(int i=0;i<m.length;i++)
				for(int j=0;j<m[i].length;j++)
					m_result[i][j] = w*m[i][j]+c;
		} catch (Exception e) {
			// TODO: handle exception
			System.out.println("出现异常，请确定输入是否正确");
			e.printStackTrace();
		}
		return m_result;
	}
	
	/**
	 * 矩阵相加
	 * @param a 输入矩阵1
	 * @param b 输入矩阵2
	 * @return 返回相加后的矩阵
	 * @throws Exception
	 */
	public static double[][] matrixAdd(double[][] a,double[][] b)
	{
		double[][] m = new double[a.length][a[0].length];
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++) 
				m[i][j]=a[i][j]+b[i][j];
		return m;
	}
	
	/**
	 * 矩阵转置
	 * @param A 输入矩阵
	 * @return 返回转置后的矩阵
	 * @throws Exception
	 */
	public static double[][] matrixTranspose(double[][] A)
	{
		double[][] AT = new double[A[0].length][A.length];
		for (int i = 0; i < AT.length; i++)
			for (int j = 0; j < AT[i].length; j++)
				AT[i][j] = A[j][i];
		return AT;
	}
	
	/**
	 * 矩阵标准相乘
	 * @param a 输入矩阵1
	 * @param b 输入矩阵2
	 * @return 返回相乘后的矩阵
	 * @throws Exception
	 */
	public static double[][] matrixProduct(double[][] a,double[][] b)
	{
		double[][] m = null;
		if(a[0].length!=b.length) m = null;
		else {
			m = new double[a.length][b[0].length];
			for (int i = 0; i < m.length; i++) {
				for (int j = 0; j < m[i].length; j++) {
					for (int k = 0;  k< b.length;k++) {
						m[i][j] +=a[i][k]*b[k][j];
					} 
				}
			}
		}
		return m;
	}
	
	/**
	 * 计算Hadamard矩阵，两矩阵元素对应相乘
	 * @param a 输入矩阵1
	 * @param b 输入矩阵2
	 * @return 返回Hadamard乘积矩阵
	 * @throws Exception
	 */
	public static double[][] hadamardProduct(double[][] a,double[][] b)
	{
		double[][] m = null;
		m = new double[a.length][a[0].length];
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				m[i][j] = a[i][j]*b[i][j];
			}
		}
		return m;
	}
	
	/**
	 * 计算行列式(m,n)坐标下的余子式
	 * 注意：m,n均从1开始,不是按照数组下标
	 * @param A 输入矩阵1
	 * @param m 横坐标
	 * @param n 纵坐标
	 * @return 返回(m,n)坐标下的余子式
	 */
	public static double[][] getCofactor(double[][] A, int m,int n)
	{
		int M= A.length;
		int N = A[0].length;
		double[][] result = new double[M-1][N-1];
		for (int i = 0; i < result.length; i++) {
			if (i<m-1) {
				for (int j = 0; j < result[i].length; j++) {
					if (j<n-1) {
						result[i][j] = A[i][j];
					}
					else result[i][j] = A[i][j+1];
				}
			}
			else{
				for (int j = 0; j < result[i].length; j++) {
					if (j<n-1) {
						result[i][j] = A[i+1][j];
					}
					else result[i][j] = A[i+1][j+1];
				}
			}
		}
		return result;
	}
	
	/**
	 * 计算行列式的值
	 * @param A 输入行列式A
	 * @return 返回行列式的值
	 */
	public static double det(double[][] A)
	{
		
		double result = 0;
		if(A.length==2) result = A[0][0]*A[1][1]-A[0][1]*A[1][0];
		else{
			double[] data = new double[A.length];
			for (int i = 0; i < data.length; i++) {
				if(i%2==0)
					data[i]=A[0][i]*det(getCofactor(A,1,i+1));
				else
					data[i]=-A[0][i]*det(getCofactor(A,1,i+1));
			}
			for (int i = 0; i < data.length; i++) {
				result+=data[i];
			}
		}
		return result;
	}
	
	/**
	 * 计算逆矩阵A-1 = A*\|A|
	 * 可逆条件：det(A)不为0，且A为方阵
	 * @param A 输入方阵A
	 * @return 返回逆矩阵，如果矩阵不可逆返回null
	 */
	public static double[][] matrixInverse(double[][] A)
	{
		double[][] result = new double[A.length][A[0].length];
		double val = det(A);
		if(val==0){
			System.out.println("行列式为0，矩阵不可逆");
			result = null;
		}
		else {
			for(int i=0; i<A.length; i++) {
	            for(int j=0; j<A[0].length; j++) {
	                if((i+j)%2 == 0) {
	                    result[i][j] = det(getCofactor(A, i+1, j+1)) /val;
	                }else {
	                    result[i][j] = -det(getCofactor(A, i+1, j+1)) /val;
	                }

	            }
	        }
			result = matrixTranspose(result);
		}
		return result;
	}
	
	/**
	 * 矩阵反转
	 * @param A 输入矩阵A
	 * @return 返回反转矩阵
	 */
	public static double[][] matrixReverse(double[][] A)
	{
		double[][] result = new double[A.length][A[0].length];
		
		for(int i=0; i<A.length; i++) {
            for(int j=0; j<A[i].length; j++) {
                result[i][j] = A[i][A[i].length-1-j];
            }
        }

		return result;
	}
	
	/**
	 * 向量反转
	 * @param A 输入向量A
	 * @return 返回反转向量
	 */
	public static double[] vectorReverse(double[] A)
	{
		double[] result = new double[A.length];
		for(int i=0; i<A.length; i++) {
            result[i] = A[A.length-1-i];
        }
		return result;
	}
	
	/**
	 * 矩阵范数
	 * @param A 输入矩阵A
	 * @param N 输入幂数N
	 * @return 返回N次范数L(N)
	 */
	public static double matrixNorm(double[][] A,int N)
	{
		double result =0;
		
		for(int i=0; i<A.length; i++) {
            for(int j=0; j<A[i].length; j++) {
                result+=Math.pow(A[i][j], N);
            }
        }
		result = Math.pow(result, 1.0/N);
		return result;
	}
	
	/**
	 * 向量范数
	 * @param A 输入向量A
	 * @param N 输入幂数N
	 * @return 返回N次范数L(N)
	 */
	public static double vectorNorm(double[] A,int N)
	{
		double result =0;
		
		for(int i=0; i<A.length; i++) {
            result+=Math.pow(A[i], N);
        }
		result = Math.pow(result, 1.0/N);
		return result;
	}
	
	/**
	 * 向量余弦
	 * @param A 输入向量A
	 * @param B 输入向量B
	 * @return 夹角余弦
	 */
	public static double vectorCosine(double[] A,double[] B)
	{
		double result =0;
		
		for(int i=0; i<A.length; i++) {
            result += A[i]*B[i];
        }
		result = result/(vectorNorm(A,2)*vectorNorm(B,2));
		return result;
	}
	
	/**
	 * 生成对角矩阵
	 * @param A 输入向量A
	 * @return 对角矩阵
	 */
	public static double[][] diagMatrix(double[] A)
	{
		double[][] result = new double[A.length][A.length];
		
		for(int i=0; i<A.length; i++) {
            result[i][i] = A[i];
        }
		return result;
	}
	
	/**
	 * 矩阵迹运算
	 * @param A
	 * @return 对角元素的和
	 */
	public static double matrixTr(double[][] A)
	{
		double result =0;
		if(A.length==A[0].length){
			for(int i=0; i<A.length; i++) {
	            result += A[i][i];
	        }
		}
		else throw new NullPointerException();
		
		return result;
	}
	public static void main(String[] args) throws Exception
	{
		double[][] a={{1,2,-1},{3,1,0},{-1,-1,-2}};
		double[] b = {0,1};
		double[] c = {1,1};
		//double[][] b={{-1,4,0},{0,1,5}};
		System.out.println( matrixTr(a));//Arrays.toString()
		System.out.println(vectorCosine(b,c));

	}
			
}
