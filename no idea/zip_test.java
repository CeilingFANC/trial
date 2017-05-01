import java.io.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;

public class zip_test {

	private static void shuffleArray(int[] array)
	{
	    int index, temp;
	    Random random = new Random();
	    int count=array.length;
	    for (int i = array.length - 1; i > 0; i--)
	    {

	        index = random.nextInt(i + 1);
	        temp = array[index];
	        array[index] = array[i];
	        array[i] = temp;
	    }
	}
    // for simply get the result by bruto force
	public static double[] naive(double[] target){
		double min=Double.POSITIVE_INFINITY;
		int res=-1;
		double temp;
		for(int i=0;i<pointlist.length;i++){
			temp=betweenpoint_2(i,target);
			if(temp<min){
				min=temp;
				res=i;
			}
		}
		return pointlist[res];
		
	}  
	
	private static double betweenpoint_2(int a, double[] b){
		double re=0;
		for (int i=0;i<num_dim;i++){
			re+=Math.pow(pointlist[a][i]-b[i], 2);
		}
		return re;
	}
	public static double[][] linear(int k,double[] target){
		if(k<1){
			return null;
		}
		int[] poten_ind=new int[k];
		double[] poten_val=new double[k];
		for (int i=0;i<k;i++){
			poten_ind[i]=-1;
			poten_val[i]=Double.MAX_VALUE;
		}
	
		double temp;
		for (int i=0;i<pointlist.length;i++){
			temp=betweenpoint_2(i,target);
			update(i,temp,poten_ind,poten_val);
			
		}
		double[][] res=new double[k][num_dim];
		for (int i=0;i<k;i++){
			res[i]=pointlist[poten_ind[i]].clone();
			
		}
		return res;
		
	}
	private static void update(int index,double dist,int[] potential_index,double[] potential_value){
		for (int i=0;i<potential_index.length;i++){
			if(dist<potential_value[i]){
				potential_value[i]=dist;
				potential_index[i]=index;
				return;
			}
		}
	}
    public static boolean check(int k,double[] target,double[] check){
    	double[][] lin;
		lin=linear(k,target);
		int res=0;
    	for (int i=0;i<k;i++){
    		res=0;
    		//System.out.println(Arrays.toString(lin[i]));
    		for (int j=0;j<num_dim;j++){
    			if(lin[i][j]==check[j]){
    				res++;
    			}
    		}
    		if(res==num_dim){
    			return true;
    		}
    	}
    	return false;
    }
    private static double count_right(){
    	int count=0;
    	for (int i=0;i<ratio_right.length;i++){
    		if(ratio_right[i]){
    			count++;
    		}
    	}
    	return (double)count/ratio_right.length;
    }
    private static double average_time1(){
    	double av_time=0;
    	for (int i=0;i<time1.length;i++){
    		av_time+=time1[i];
    	}
    	return av_time/time1.length;
    }
    private static double average_time2(){
    	double av_time=0;
    	for (int i=0;i<time2.length;i++){
    		av_time+=time2[i];
    	}
    	return av_time/time2.length;
    }    
    
    
    
    
    
    

	public static int selectKth(int[] arr,int dim,int from,int to) {
		if(from==to){
			 return from;
		}
		int tmp;
		//System.out.print(to);

		int k=(from+to)/2;
		 // if from == to we reached the kth element
		while (from < to) {
			int r = from;
			int w = to;
			double mid = pointlist[arr[(r + w) / 2]][dim];
			//System.out.print(mid);
	 
		  // stop if the reader and writer meets
			while (r < w) {
				if (pointlist[arr[r]][dim] >= mid) { // put the large values at the end
					tmp = arr[w];
					arr[w] = arr[r];
					arr[r] = tmp;
					w--;
					//System.out.println("no"+w);
				} else { // the value is smaller than the pivot, skip
					r++;
					//System.out.println("yes"+r);

				}
			}
		  // if we stepped up (r++) we need to step one down
			if (pointlist[arr[r]][dim] > mid){
				r--;
			}
		  // the r pointer is on the end of the first k elements
			if (k <= r) {
				to = r;
	
			} else {
				from = r + 1;
			}

		 }
		 
		 return k;
	}
    private static void selectionSort(Kdtree tr,int[] inse,int dim,int from,int to){ 
    	if (to-from>5){
    		return;
    	}
        for (int i = from; i < to+1; i++)  
        {  
            int index = i;  
            for (int j = i + 1; j < to+1; j++){  
                if (pointlist[inse[j]][dim] < pointlist[inse[index]][dim] ){  
                    index = j;//searching for lowest index  
                }  
            }  
            int smallerNumber = inse[index];   
            inse[index] = inse[i];  
            inse[i] = smallerNumber;  
        }
		Point2D ref1=new Point2D(pointlist[inse[(to+from)/2]][0],pointlist[inse[(to+from)/2]][1]);
		ref1.link=inse[from];
        tr.insert(ref1);
        for (int i=from;i<to+1;i++){
        	if(i!=(to+from)/2){
	    		Point2D ref2=new Point2D(pointlist[inse[i]][0],pointlist[inse[i]][1]);
	    		ref2.link=inse[i];
        		tr.insert(ref2);
        	}
        }
    }  
	private static void order(Kdtree tr,int[] inse, int dim,int from, int to){

		if(to-from<5){
			selectionSort(tr,inse,dim,from,to);
//			for (int i=from;i<to+1;i++){
//	    		Point2D ref=new Point2D(pointlist[inse[i]][0],pointlist[inse[i]][1]);
//	    		ref.link=inse[i];
//				tr.insert(ref);
//			}
		}else{
			int split_1=selectKth(inse,rotate[dim],from,to);
			if(to!=-1){
				//System.out.println(cooo+"^"+split_1+"which one!!!"+from+"--"+to+"__"+remember_last);
			}
    		Point2D ref=new Point2D(pointlist[inse[split_1]][0],pointlist[inse[split_1]][1]);
    		ref.link=inse[split_1];
			tr.insert(ref);
			if(split_1-1==from){
	    		Point2D ref1=new Point2D(pointlist[inse[from]][0],pointlist[inse[from]][1]);
	    		ref1.link=inse[from];
				tr.insert(ref1);
			}else{
				order(tr,inse,(dim+1)%use_dim,from,split_1-1);
			}
			if(split_1+1==to){
	    		Point2D ref2=new Point2D(pointlist[inse[to]][0],pointlist[inse[to]][1]);
	    		ref2.link=inse[to];
				tr.insert(ref2);
			}else{
				order(tr,inse,(dim+1)%use_dim,split_1+1,to);
			}
		}
	}
    
    
    
	private static void mean(){
		Double sum=0.0;
		for (int i=0;i<num_dim;i++){
			sum=0.0;
			for (int j=0;j<pointlist.length;j++){
				sum+=pointlist[j][i];
			}
			//System.out.println(sum);
			means[i]=sum/pointlist.length;
		}
	}
	private static void variance(){
		Double sum=0.0;
		for (int i=0;i<num_dim;i++){
			sum=0.0;
			System.out.println(sum);

			double mea=means[i];
			for (int j=0;j<pointlist.length;j++){
				sum+=Math.pow(pointlist[j][i]-mea, 2);
				//if(i==1&&Math.pow(pointlist[j][i]-mea, 2)<0)
				//System.out.println(Math.pow(pointlist[j][i]-mea, 2));
			}
			//System.out.println(sum*1000);

			variances[i]=sum;
		}		
	}
    
    
    
    
    
    
    
    private double cha_x(double x){
    	return (x+15)/100;
    }
    private double cha_y(double y){
    	return (y+180);
    }
	private static int use_dim=2;
	private static int num_dim=2;
	private static String[] place= new String[215956/5];
	private static int[] zip=new int[215956/5];
	private static double[] latitude=new double[215956/5];
	private static double[] longtitude=new double[215956/5];
	private static double[][] pointlist;
	private static int[] rotate=new int[] {0,1};
	
	private static boolean[] ratio_right;
	private static double[] time1;
	private static double[] time2;
	
	private static double[] means=new double[] {0,0};
	private static double[] variances=new double[] {0,0};

	public static void main(String args[]) throws FileNotFoundException{

        Scanner scanner = new Scanner(new File("zipcode.csv"));
        scanner.useDelimiter(",");
        String[][] aa=new String[43195][7];
        int incre=0;
        // intialize the array for dataset
        pointlist=new double[215956/5][2];
        while(scanner.hasNext()){

        	if(incre%5==1){

        		zip[incre/5]=Integer.parseInt(scanner.next());
        	}else if(incre%5==2){
        		place[incre/5]=scanner.next();
        	}else if(incre%5==3){
        		latitude[incre/5]=Double.parseDouble(scanner.next());
        		pointlist[incre/5][0]=latitude[incre/5];
        	}else if(incre%5==4){
        		longtitude[incre/5]=Double.parseDouble(scanner.next());
        		pointlist[incre/5][1]=longtitude[incre/5];
        	} else{
        		String cc=scanner.next();
        	}
        	incre+=1;
        }
        scanner.close();
        //System.out.println(incre);
        //System.out.println(place[place.length-1]);
        
        
        // -10 to 80 for pointlist[0] latitude
        // - 55 to -180
        
        int number=500000;
        pointlist=new double[number][2];
    	Random k = new Random(); 
        for (int i=0;i<number;i++){
        	pointlist[i][0]= (-10)+90*k.nextDouble();
        	pointlist[i][1]= (-180)+125*k.nextDouble();
        }
 
        mean();
        variance();
        System.out.println(Arrays.toString(means));

        System.out.println(Arrays.toString(variances));
        
        
        
    	Kdtree million = new Kdtree();
    	int times=20000;
		time1=new double[times];
		time2=new double[times];
		ratio_right=new boolean[times];
        double[][] test=new double[times][2];
    	Random r = new Random(); 
        int[] rand=new int[pointlist.length];
        for(int i=0;i<pointlist.length;i++){
        	rand[i]=i;
        }
        shuffleArray(rand);
        //System.out.println(Arrays.toString(rand));
        order(million,rand,0,0,pointlist.length-1);
        
        
//    	for (int i=0; i<pointlist.length;i++){
//    		Point2D ref=new Point2D(pointlist[rand[i]][0],pointlist[rand[i]][1]);
//    		ref.link=rand[i];
//    		million.insert(ref);
//    	}
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	

      
        for (int i=0;i<test.length;i++){
        	test[i][0]= (-10)+90*r.nextDouble();
        	test[i][1]= (-180)+125*r.nextDouble();
        }   	

    	for (int i=0;i<test.length;i++){
    		Point2D ty=new Point2D(test[i][0],test[i][1]);
    		Point2D res;
        	long startTime = System.nanoTime();
        	res=million.nearest(ty);
        	time2[i] = (System.nanoTime() - startTime) / 1000000.0;

        	ratio_right[i]=check(5,new double[] {test[i][0],test[i][1]},new double[] {res.x(),res.y()});
    	}
    	
    	for (int i=0;i<test.length;i++){
        	long startTime = System.nanoTime();
    		//linear(1,test[i]);
        	naive(test[i]);
        	time1[i] = (System.nanoTime() - startTime) / 1000000.0;
    	}
    	double av_time1=average_time1();
    	double av_time2=average_time2();
    	double rat=count_right();
   	
    	System.out.println(av_time1);// linear
    	System.out.println(av_time2);// kd
    	System.out.println(rat);
    	//System.out.println(Arrays.toString(ratio_right));

        
        

	}
	
}





























