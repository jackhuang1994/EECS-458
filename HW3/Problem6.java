// Program: Problem6.java
// Author: Mingan Huang 
// Network ID: mxh805
// A program calculate viterbi path posterior probability 

import java.util.List;
import java.util.ArrayList;
import java.lang.*; 
import java.io.File;
import java.util.Scanner;
import java.io.FileNotFoundException;

public class Problem6{  

  static int[] states = new int[]{0, 1}; // 0 for Fair, 1 for Biased
  static double[] start_probability = new double[]{0.5, 0.5}; // 0.5 to choose a fair coin, 0.5 to choose a biased coin
  static double[][] transititon_probability = new double[][]{
    {0.9, 0.1}, // Fair to Fair, Fair to Biased
    {0.1, 0.9}, // Biased to Fair, Biased to Biased
  };
  static double[][] emission_probability = new double[][]{
    {0.5, 0.5}, // Fair: Head, Fair: Tail
    {0.75, 0.25}, // Biased: Head, Biased: Tail
  };
  
  // obs: the observation we get from input file
  // states: hidden states
  // start: start probabilities
  // trans: transtion probabilities
  // emit: emmision probabilities 
  public static int[] viterbi(int[] obs, int[] states, double[] start, double[][] trans, double[][] emit)
  {
    double[][] matrix = new double[obs.length][states.length];
    int[][] path = new int[states.length][obs.length];

    for (int x : states)
    {
      matrix[0][x] = start[x] * emit[x][obs[0]];
      path[x][0] = x;
    }

    for (int i = 1; i < obs.length; i++)
    {
      int[][] newpath = new int[states.length][obs.length];

      for (int x : states)
      {
        double prob = -1;
        int state;
        for (int y : states)
        {
          double nprob = matrix[i - 1][y] * trans[y][x] * emit[y][obs[i]];
          if (nprob > prob)
          {
            prob = nprob;
            state = y;
            // Record the maximum probability
            matrix[i][x] = prob;
            // Record the path
            System.arraycopy(path[state], 0, newpath[x], 0, i);
            newpath[x][i] = x;
          }
        }
      }

      path = newpath;
    }

    double prob = -1;
    int state = 0;
    for (int x : states)
    {
      if (matrix[obs.length - 1][x] > prob)
      {
        prob = matrix[obs.length - 1][x];
        state = x;
      }
    }
    
    System.out.println();
    System.out.println("Viterbi matrix: "); 
    for(int i = 0; i < matrix.length; i++) // print out the viterbi matrix
    {
      for (int j = 0; j < matrix[i].length; j++)
        System.out.print(matrix[i][j]+ "\t");
      System.out.println();
    }
    
    return path[state];
  }
  
  // same parameters as viterbi
  // calculate and print out the forward matrix
  public static double[][] forward_matrix(int[] obs, int[] states, double[] start, double[][] trans, double[][] emit)
  {
    double[][] matrix = new double[obs.length][states.length];
    
    for (int x : states)
    {
      matrix[0][x] = start[x] * emit[x][obs[0]];
    }
    
    for (int i = 1; i < obs.length; i++)
    {
      for (int x : states)
      {
        for (int y : states)
        {
          matrix[i][x] += emit[x][obs[i]] * trans[y][x] * matrix[i-1][y];
          // Based on the definition of forward algorithm
        }
      }
    }
    
    System.out.println();
    System.out.println("forward matrix: ");
    printMatrix(matrix);
    
    return matrix;
  }
  
  // same parameters as viterbi
  // calculate and print out the backward matrix
  public static double[][] backward_matrix(int[] obs, int[] states, double[] start, double[][] trans, double[][] emit)
  {
    double[][] matrix = new double[obs.length][states.length];
    
    for (int x : states)
    {
      matrix[obs.length-1][x] = 1;
    }
    
     for (int i = 1; i < obs.length; i++)
    {
      for (int x : states)
      {
        for (int y : states)
        {
          matrix[obs.length - i - 1][y] += emit[x][obs[obs.length - i - 1]] * trans[x][y] * matrix[obs.length - i][y];
          // Based on the definition of backward algorithm
        }
      }
    }
    
    System.out.println("backward matrix: ");
    printMatrix(matrix);
    
    return matrix;
  }
  
  // Using the forward and backward matrixes to calculate posterior probability
  public static double posterior_probability(int pos, int state, double[][] forward_matrix, double[][] backward_matrix)
  {
    double result = 0;
    int length = forward_matrix.length;
    
    double forward_state = forward_matrix[pos][state];
    double backward_state = backward_matrix[pos][state];
    result = forward_state * backward_state / (forward_matrix[length-1][0] + forward_matrix[length-1][1]);
    // Using information in forward matrix and backward matrix to calculate posterior probability
    
    return result;
  }
  
  // Prints out the whole matrix
  public static void printMatrix(double[][] matrix) {
    for ( int i = 0; i < matrix.length; i++ ) {
      for ( int j = 0; j < matrix[i].length; j++ ) 
        System.out.print(matrix[i][j] + "\t");
      System.out.println();
    }
    System.out.println();
  }
  
  public static void main(String[] args) throws FileNotFoundException {
    String obs = "";
    int position = 0;    
    
    System.out.println("Please input the name of the file: (do not include .txt at the end! e.g. input1)");
    Scanner input = new Scanner(System.in);
    String file = (input.nextLine() + ".txt");
      
    try {
      Scanner fileInput = new Scanner(new File(file));
      while(fileInput.hasNextLine()) {
        String newLine = fileInput.nextLine();
        if (newLine.toLowerCase().startsWith("h") || newLine.toLowerCase().startsWith("t")) {
          obs += newLine; // Store the sequence in obs
        }
        else 
        {
          String[] splitStr = newLine.split(" ");
          position = Integer.parseInt(splitStr[0]); // store the position for posterior probability in position
        }
      }
    } catch (FileNotFoundException e) {
      System.err.println("File not found");
    }

    int[] observations = new int[obs.length()];
    for(int i = 0; i < obs.length(); i++) // convert string obs to an array of integer
    {
      Character obj = new Character('H'); 
      if (obs.charAt(i) == obj)
        observations[i] = 0;
      else
        observations[i] = 1;
    }

    int[] result = viterbi(observations, states, start_probability, transititon_probability, emission_probability);
    
    System.out.println();
    System.out.println("Viterbi path: ");
    for (int r : result) // print out the viterbi path
    {
      if (r == 0)
        System.out.print("Fair ");
      else
        System.out.print("Biased ");
    }
    System.out.println();
    
    double[][] forward = forward_matrix(observations, states, start_probability, transititon_probability, emission_probability);
    double[][] backward = backward_matrix(observations, states, start_probability, transititon_probability, emission_probability);
    double posterior = posterior_probability(position, 1, forward, backward); 
    // Posterior probability of position 7 generated by a biased coin, since 1 represent states biased
    System.out.println("The posterior probability of position " + position + " is generated by a biased coin is " + posterior);
  }
}
  