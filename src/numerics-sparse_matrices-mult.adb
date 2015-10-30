separate (Numerics.Sparse_Matrices)

function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is
   use Ada.Containers, Ada.Text_IO;
   A  : Sparse_Matrix renames Left;
   B  : Sparse_Matrix renames Right;
   C  : Sparse_Matrix;
   X  : Real_Array (1 .. A.N_Row) := (others => 0.0);
   W  : Int_Array (1 .. A.N_Row)  := (others => 0);
   Nz : Pos := 1;
   
   procedure Scatter (A	   : in     Sparse_Matrix;
		      J	   : in     Int;
		      β	   : in     Real;
		      W	   : in out Int_Array;
		      X	   : in out Real_Array;
		      Mark : in     Int;
		      C	   : in out Sparse_Matrix;
		      Nz   : in out Int) is
      I : Int;
   begin
      for P in A.P (J) .. A.P (J + 1) - 1 loop
	 I := A.I (P);
	 if W (I) < Mark then
	    W (I) := Mark;
	    C.I.Append (I);
	    Nz := Nz + 1;
	    X (I) := β * A.X (P);
	 else
	    X (I) := X (I) + β * A.X (P);
	 end if;
      end loop;
   end Scatter;
   
begin
   C.Format := CSC; C.N_Row := A.N_Row; C.N_Col := B.N_Col;
   C.P.Reserve_Capacity (Count_Type (B.N_Col) + 1);
   
   for J in 1 .. B.N_Col loop
      C.P.Append (Nz);
      for K in B.P (J) .. B.P (J + 1) - 1 loop
	 Scatter (A, B.I (K), B.X (K), W, X, J, C, Nz);
      end loop;
      
      for P in C.P (J) .. Nz - 1 loop
	 C.X.Append (X (C.I (P)));
      end loop;
   end loop;
   C.P.Append (Nz);
   
   C.Convert; 
   C.Convert;
   return C;
end Mult;
