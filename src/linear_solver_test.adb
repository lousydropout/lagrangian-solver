with Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
use  Numerics, Numerics.Sparse_Matrices, Numerics.Sparse_Matrices.CSparse;
with Ada.Text_IO; use Ada.Text_IO;

procedure Linear_Solver_Test is
   use Real_IO, Int_IO, Real_Functions;
   
   Mat : Sparse_Matrix;
   LU  : LU_Type;
   X   : Real_Vector;
   B   : Real_Vector;
   Dir : String      := "matrices/sparse-triplet/zero-based/";
   Res : Real;

begin
   Put ("Read Matrix . . .");
   
   --  Mat := Read_Sparse_Triplet (Dir & "a5by5_st.txt");     -- 611B
   --  Mat := Read_Sparse_Triplet (Dir & "bcsstk01_st.txt");  -- 4.9K
   Mat := Read_Sparse_Triplet (Dir & "bcsstk16_st.txt");  -- 3.7M
   --  Mat := Read_Sparse_Triplet (Dir & "fs_183_1_st.txt");  -- 24K
   --  Mat := Read_Sparse_Triplet (Dir & "kershaw_st.txt");   -- 564
   --  Mat := Read_Sparse_Triplet (Dir & "west0067_st.txt");  -- 3.9K
   --  Mat := Read_Sparse_Triplet (Dir & "t1_st.txt");        -- 80
   Put_Line ("finished");
   

   ----- Print matrix' info --------------
   Put ("Size of matrix: "); 
   Put (Mat.N_Row, 0); Put (" x "); Put (Mat.N_Col, 0); New_Line;
   Put ("Number of entries: "); Put (Mat.Number_Of_Elements, 0); New_Line;
   
   ----- Set size of vectors X and B ----
   Set_Length (B, Mat.N_Col); Set_Length (X, Mat.N_Col);

   ----- Begin LU Decomposition ---------
   Put ("LU Decomposition . . .");
   LU  := LU_Decomposition (Mat);
   Put_Line ("finished");
   
   ------ Begin tests ------------------------
   Put_Line ("Begin testing . . .");
   for K in 1 .. Int (10) loop
      Put ("Trial "); Put (K, Width => 2); Put (": "); 
      for I of B loop I := (10.0 * Rand) ** 10 * Sin (10.0 * Rand); end loop;
      X   := Solve (LU, B);
      Res := Norm (B - Mat * X) / (if Norm (X) > 1.0 then Norm (X) else 1.0);
      Put ("    Norm (Res)  =  "); Put (Res, Fore => 1, Aft => 1, Exp => 3);
      if Res > 1.0e-10 then Put ("  ***"); end if;
      New_Line;
   end loop;
   Put_Line ("tests completed");
   
   ---- Free memory allocated by LU : LU_Type -------
   Free (LU);
   
end Linear_Solver_Test;
