with Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
with Ada.Text_IO; use Ada.Text_IO;
with Ada.Containers; use Ada.Containers;
procedure Linear_Solver_Test is
   use Real_IO, Int_IO, Real_Functions;
   
   Mat : Sparse_Matrix;
   LU  : LU_Type;
   SX  : Sparse_Vector;
   SB  : Sparse_Vector;
   SY  : Sparse_Vector;
   Y   : Sparse_Vector;
   Dir : String      := "matrices/sparse-triplet/zero-based/";
   Res : Real;

begin
   Put ("Read Matrix . . .");
   
   ------ Rectangular Matrices -----------------
   --  Mat := Read_Sparse_Triplet (Dir & "ash219_st.txt");       -- 3.6K
   --  Mat := Transpose (Mat) * Mat;
   
   ------ Square Matrices ----------------------
   --  Mat := Read_Sparse_Triplet (Dir & "a5by5_st.txt");     -- 611
   --  Mat := Read_Sparse_Triplet (Dir & "bcsstk01_st.txt");  -- 4.9K
   Mat := Read_Sparse_Triplet (Dir & "bcsstk16_st.txt");  -- 3.7M
   --  Mat := Read_Sparse_Triplet (Dir & "fs_183_1_st.txt");  -- 24K
   --  Mat := Read_Sparse_Triplet (Dir & "kershaw_st.txt");   -- 564
   --  Mat := Read_Sparse_Triplet (Dir & "t1_st.txt");        -- 80
   --  Mat := Read_Sparse_Triplet (Dir & "west0067_st.txt");  -- 3.9K
   
   Put_Line ("finished");

   
   ----- Print matrix' info --------------
   Put ("Size of matrix: "); 
   Put (Mat.N_Row, 0); Put (" x "); Put (Mat.N_Col, 0); New_Line;
   Put ("Number of entries: "); Put (Mat.Number_Of_Elements, 0); New_Line;
   
   
   ----- Set size of vectors X and B ----
   Set_Length (SB, Mat.N_Col); Set_Length (SX, Mat.N_Col); 
   
   
   ----- Begin LU Decomposition ---------
   Put ("LU Decomposition . . .");
   LU  := LU_Decomposition (Mat);
   Put_Line ("finished");
   
   
   ------ Begin tests ------------------------
   Put_Line ("Begin testing . . .");
   for K in 1 .. Int (10) loop
      Put ("Trial "); Put (K, Width => 2); Put (": "); 
      
      for I in 1 .. Mat.N_Col loop
	 SB.Set (I, (10.0 * Rand) ** 10 * Sin (10.0 * Rand));
      end loop;
      
      --  SX  := Solve (LU, SB);
      --  SY  := SB - Mat * SX;
      Y   := SB - Mat * Solve (Mat, SB);
      --  Res := Norm (SY - Y);
      Res := Norm (Y) / Real'Max (1.0, Norm (SB));
      Put ("    Norm (Res)  =  "); Put (Res, Fore => 1, Aft => 1, Exp => 3);
      
      Put_Line (if Res > 1.0e-10 then "  ***" else "");
   end loop;
   Put_Line ("tests completed");
   
   ---- Free memory allocated by LU : LU_Type -------
   Free (LU);
   
end Linear_Solver_Test;
