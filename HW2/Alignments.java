// Program: Alignments.java
// Author: Mingan Huang 
// Network ID: mxh805
// A function to do both global and local alignment for two sequences 
// Modified some functions from internet

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.util.Scanner;
import java.io.FileNotFoundException;

public class Alignments {  
  public static int globalAlignment(String sequence1, String sequence2, int match, int mismatch, int indel) {
    int m = sequence1.length();
    int n = sequence2.length();
 
    int[][] scores = new int[m+1][n+1];  // create socres table for optimal score
    int[][] backtracking = new int[m+1][n+1];  // create backtracking table for back trackpath
    
    // uses 3 bit binary (represented as an int) to encode directions
    // up, left, diagonal in that order (e.g. 011 represents to left and diagonal)
 
    // initialize score matrix
    for ( int i = 0; i <= m; i++ )
      scores[i][0] = indel * i;
    for ( int j = 0; j <= n; j++ )
      scores[0][j] = indel * j;

    // initialize backtrack matrix
    backtracking[0][0] = 0;
    for ( int i = 1; i <= m; i++ ) 
      backtracking[i][0] = 4;
    for ( int j = 1; j <= n; j++ ) 
      backtracking[0][j] = 2;
 
 
    // fill in subproblem matrices
    for ( int i = 1; i <= m; i++ ) {
      for ( int j = 1; j <= n; j++ ) {
        int sub = scores[i-1][j-1] + sub(sequence1.charAt(i-1), sequence2.charAt(j-1), match, mismatch); // score for substitution
        int gap2 = scores[i-1][j] + indel; // score for gap in s2
        int gap1 = scores[i][j-1] + indel; // score for gap in s1

        // fill in solution
        int maximum = max(sub, gap2, gap1);
        scores[i][j] = maximum;

        // fill in backtrack data
        if ( maximum == sub ) // substitution gave the best optimal subsolution
          backtracking[i][j] += 1; // point diagonal
        if ( maximum == gap2 ) // inserting gap in s2 gave the best optimal subsolution
          backtracking[i][j] += 4; // point up
        if ( maximum == gap1 ) // inserting gap in s1 gave the best optimal subsolution
          backtracking[i][j] += 2; // point left
      }

    }

 
    printMatrix(scores);
    printMatrix(backtracking);
    System.out.println("Optimal score: " + scores[m][n]);
    System.out.println("Number of solutions: " + alignmentsCount(false, scores, backtracking, m, n));
    System.out.println("Actual Aligenments: ");
    alignmentsPrint(false, sequence1, sequence2, scores, backtracking, m, n);
 
    return scores[m][n];
  }
  
  public static int localAlignment(String sequence1, String sequence2, int match, int mismatch, int indel) {
    // length of input strings
    int m = sequence1.length();
    int n = sequence2.length();

    int[][] scores = new int[m+1][n+1];    // create score solution matrix
    int[][] backtracking = new int[m+1][n+1]; // create backtracking matrix
    
    // uses 3 bit binary (represented as an int) to encode directions
    // up, left, diagonal in that order (e.g. 011 represents to left and diagonal)
 
    // initialize score matrix (same as Needleman-Wunsch)
    for ( int i = 0; i <= m; i++ ) 
      scores[i][0] = indel * i;
    for ( int j = 0; j <= n; j++ ) 
      scores[0][j] = indel * j;
    
    // initialize backtrack matrix (different from Needleman-Wunsch)
    backtracking[0][0] = 0;
    for ( int i = 1; i <= m; i++ ) 
      backtracking[i][0] = 0;
    for ( int j = 1; j <= n; j++ ) 
      backtracking[0][j] = 0;

    // fill in subproblem matrices
    for ( int i = 1; i <= m; i++ ) {
      for ( int j = 1; j <= n; j++ ) {
        int sub = scores[i-1][j-1] + sub(sequence1.charAt(i-1), sequence2.charAt(j-1), match, mismatch); // score for substitution
        int gap2 = scores[i-1][j] + indel; // score for gap in s2
        int gap1 = scores[i][j-1] + indel; // score for gap in s1
  
        // fill in solution
        int maximum = max(sub, gap2, gap1);
        // best score is 0 if we got a negative value
        if ( maximum < 0 ) 
          scores[i][j] = 0;
        else 
          scores[i][j] = maximum;

        // fill in backtrack data
        if ( scores[i][j] == sub ) // substitution gave the best optimal subsolution
          backtracking[i][j] += 1; // point diagonal
        if ( scores[i][j] == gap2 ) // inserting gap in s2 gave the best optimal subsolution
          backtracking[i][j] += 4; // point up
        if ( scores[i][j] == gap1 ) // inserting gap in s1 gave the best optimal subsolution
          backtracking[i][j] += 2; // point left
      }
    }

    // find optimal score in matrix
    int optimalScore = scores[0][0];
    
    // the indices (i, j) of the optimal solutions (there can be more than one)
    for ( int i = 0; i <= m; i++ ) {
      for ( int j = 0; j <= n; j++ ) {
        if ( scores[i][j] > optimalScore ) 
          optimalScore = scores[i][j];
      }
    }

    // find the indices where the max scores are
    List<int[]> maxScore = new ArrayList<int[]>();
    for ( int i = 0; i <= m; i++ ) {
      for ( int j = 0; j <= n; j++ ) {
        if ( scores[i][j] == optimalScore )
          maxScore.add(new int[] {i, j});
      }
    }
 
    // count the number of optimal alignments
    int numOptimalAlignments = 0;
    for ( int k = 0; k < maxScore.size(); k++ ) {
      int i = maxScore.get(k)[0];
      int j = maxScore.get(k)[1];
      numOptimalAlignments += alignmentsCount(true, scores, backtracking, i, j);
    }
 
    printMatrix(scores);
    printMatrix(backtracking);
    
    System.out.println("Optimal score: " + optimalScore);
    System.out.println("Number of solutions: " + numOptimalAlignments);
    System.out.println("Actual Aligenments: ");
    
    // print all the optimal alignments
    for ( int k = 0; k < maxScore.size(); k++ ) {
      int i = maxScore.get(k)[0];
      int j = maxScore.get(k)[1];
      alignmentsPrint(true, sequence1, sequence2, scores, backtracking, i, j);
    }

    return optimalScore;
  }
  
  // returns the maximum value of three inputs
  private static int max(int x, int y , int z) {
    if ( x >= y && x >= z )
      return x;
    else if ( y >= x && y >= z )
      return y;
    else
      return z;
  }

  // Takes four inputs: chars x,y and ints match, mismatch
  // Returns the cost of match if chars match, cost of mismatch if chars mismatch
  public static int sub(char x, char y, int match, int mismatch) {
    if ( x == y )
      return match;
    else
      return mismatch;
  }
  
  // first input indicates it is getting a global or local alignment
  // Takes the sequence Strings and matrix representing backtracking data
  // Also takes a starting index to backtrack from
  // Prints out all the paths to the optimal solution Recursively
  public static List<String[]> getAlignments(boolean l, String sequence1, String sequence2, int[][] scores, int[][] backtrack, int i, int j) {
    if ( i == 0 & j == 0 ) { // base case
      List<String[]> newList = new ArrayList<String[]>();
      newList.add(new String[] {"", ""});
      return newList;
    }

    List<String[]> newAlignments = new ArrayList<String[]>();
    // if (i, j) points diagonal
    if ( backtrack[i][j] %  2 == 1 ) {
      // get all the optimal alignment sequences so far to index (i-1, j-1)
      List<String[]> setOfAlignments = getAlignments(l, sequence1, sequence2, scores, backtrack, i-1, j-1);
      for ( int k = 0; k < setOfAlignments.size(); k++ ) { // for all alignments in setsAlignments
        String[] alignment = setOfAlignments.get(k);
        alignment[0] = alignment[0] + sequence1.charAt(i-1);
        alignment[1] = alignment[1] + sequence2.charAt(j-1);
        newAlignments.add(alignment);
      }
    }

    // if (i, j) points left
    if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 ) {
      // get all the optimal alignment sequences so far to index (i, j-1)
      List<String[]> setsAlignments = getAlignments(l, sequence1, sequence2, scores, backtrack, i, j-1);
      for ( int k = 0; k < setsAlignments.size(); k++ ) {
        String[] alignment = setsAlignments.get(k);
        alignment[0] = alignment[0] + "-";
        alignment[1] = alignment[1] + sequence2.charAt(j-1);
        newAlignments.add(alignment);
      }
    }

    // if ( i, j) points up
    if ( backtrack[i][j] >= 4 ) {
      // get all the optimal alignment sequences so far to index (i-1, j)
      List<String[]> setsAlignments = getAlignments(l, sequence1, sequence2, scores, backtrack, i-1, j);
      for ( int k = 0; k < setsAlignments.size(); k++ ) {
        String[] alignment = setsAlignments.get(k);
        alignment[0] = alignment[0] + sequence1.charAt(i-1);
        alignment[1] = alignment[1] + "-";
        newAlignments.add(alignment);
      }
    }

    // if (i, j) can also be a new alignment
    if ( l && scores[i][j] == 0 ) {
      newAlignments.add(new String[] {"", ""});
    }
    
    return newAlignments;
  }
  
  // first input indicates it is getting a global or local alignment
  // recursively Counts the number of optimal alignments
  public static int alignmentsCount(boolean l, int[][] scores, int[][] backtrack, int i, int j) {
    // base case
    if ( backtrack[i][j] == 0 ) // end of sequence
      return 1;
 
    int countAlignments = 0;
    // if (i, j) points diagonal
    if ( backtrack[i][j] % 2 == 1 )
      // count alignment sequences in the diagonal direction
     countAlignments += alignmentsCount(l, scores, backtrack, i-1, j-1);
    // if (i, j) points left
    if ( backtrack[i][j] >= 2 && backtrack[i][j] / 2 != 2 )
      // count alignment sequences in the left direction
      countAlignments += alignmentsCount(l, scores, backtrack, i, j-1);
    // if (i, j) points up
    if ( backtrack[i][j] >= 4 )
      // count alignment sequences in the up direction
     countAlignments += alignmentsCount(l, scores, backtrack, i-1, j);

    // if (i, j) can also be the start of a new sequence
    if ( l && scores[i][j] == 0 )
      countAlignments += 1;

    return countAlignments;
  }
  
  // first input indicates it is getting a global or local alignment
  // prints all alignments one by one
   public static void alignmentsPrint(boolean l, String s1, String s2, int[][] scores, int[][] backtrack, int i, int j) {
    List<String[]> alignments = getAlignments(l, s1, s2, scores, backtrack, i, j);
    for ( int m = 0; m < alignments.size(); m++ ) {
      System.out.println(alignments.get(m)[0]);
      System.out.println(alignments.get(m)[1]);
      System.out.println();
    }
    System.out.println();
  }
    
  // Prints out the whole matrix
  public static void printMatrix(int[][] matrix) {
    for ( int i = 0; i < matrix.length; i++ ) {
      for ( int j = 0; j < matrix[i].length; j++ ) 
        System.out.print(matrix[i][j] + "\t");
      System.out.println();
    }
    System.out.println();
  }
    
  public static void main(String[] args) throws FileNotFoundException{
    String alignment = "";
    int match = 0;
    int mismatch = 0;
    int indel = 0;
    String sequence1 = "";
    String sequence2 = "";
    List<String> sequences = new ArrayList<String> ();
    
    System.out.println("Please input the name of the file: (including .txt at the end)");
    Scanner input = new Scanner(System.in);
    String file = input.nextLine();
      
    try {
      Scanner fileInput = new Scanner(new File(file));
      while(fileInput.hasNextLine()) {
        String newLine = fileInput.nextLine();
        if (newLine.toLowerCase().startsWith("g") || newLine.toLowerCase().startsWith("l")) {
          alignment += newLine; // store "g" in alignment if it is a global alignment, "l" if it is a local alignment 
        }
        else if (newLine.startsWith("1")) {
          String[] splitStr = newLine.split(" ");
          match = Integer.parseInt(splitStr[0]); // read score for match
          mismatch = Integer.parseInt(splitStr[1]); // read score for mismatch
          indel = Integer.parseInt(splitStr[2]); // read socre for indel
        }
        else 
          sequences.add(newLine); // add last to sequences to a list of string
      }
      sequence1 = sequences.get(0); // get the first sequence
      sequence2 = sequences.get(1); // get the second sequence
    } catch (FileNotFoundException e) {
      System.err.println("File not found");
    }
    
    if (alignment.compareTo("g") == 0) {
      System.out.println("Doing global alignment...");
      System.out.println("Match score: " + match);
      System.out.println("Mismatch score: " + mismatch);
      System.out.println("Indel score: " + indel);
      globalAlignment(sequence1, sequence2, match, mismatch, indel);
    }
    else {
      System.out.println("Doing local alignment...");
      System.out.println("Match score: " + match);
      System.out.println("Mismatch score: " + mismatch);
      System.out.println("Indel score: " + indel);
      localAlignment(sequence1, sequence2, match, mismatch, indel);
    }
  }
}
