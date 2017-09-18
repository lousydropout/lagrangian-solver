with Numerics, Ada.Text_IO;
use  Numerics;

package body Numerics.Dense_Matrices is
   
   function Outer (X, Y : in Real_Vector) return Real_Matrix is
      Result : Real_Matrix (1 .. X'Length, 1 .. Y'Length);
   begin
      for I in Result'Range (1) loop
	 for J in Result'Range (2) loop
	    Result (I, J) := X (I) * Y (J);
	 end loop;
      end loop;
      return Result;
   end Outer;
   
   function "-" (A : in Real_Matrix) return Real_Matrix is
      Result : Real_Matrix := A;
   begin
      for X of Result loop X := -X; end loop;
      return Result;
   end "-";
   
   function "*" (X : in Real;
		 A : in Real_Matrix) return Real_Matrix is
      Result : Real_Matrix := A;
   begin
      for Y of Result loop Y := X * Y; end loop;
      return Result;
   end "*";

   function "+" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix is
      C : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2))
	:= (others => (others => 0.0));
      Coff_1 : constant Integer := 1 - A'First (1);
      Coff_2 : constant Integer := 1 - A'First (2);
      Boff_1 : constant Integer := B'First (1) - A'First (1);
      Boff_2 : constant Integer := B'First (2) - A'First (2);
   begin
      for I in A'Range (1) loop
	 for J in A'Range (2) loop
	    C (I + Coff_1, J + Coff_2) := A (I, J) + B (I + Boff_1, J + Boff_2);
	 end loop;
      end loop;
      return C;
   end "+";
      
   function "-" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix is
      C : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2))
	:= (others => (others => 0.0));
      Coff_1 : constant Integer := 1 - A'First (1);
      Coff_2 : constant Integer := 1 - A'First (2);
      Boff_1 : constant Integer := B'First (1) - A'First (1);
      Boff_2 : constant Integer := B'First (2) - A'First (2);
   begin
      for I in A'Range (1) loop
	 for J in A'Range (2) loop
	    C (I + Coff_1, J + Coff_2) := A (I, J) - B (I + Boff_1, J + Boff_2);
	 end loop;
      end loop;
      return C;
   end "-";
      
   function "*" (A : in Real_Matrix;
		 X : in Real_Vector) return Real_Vector is
      B : Real_Vector (1 .. A'Length (1)) := (others => 0.0);
      Offset_1 : constant Integer := A'First (1) - 1;
      Offset_2 : constant Integer := A'First (2) - X'First;
   begin
      for I in B'Range loop
	 for J in X'Range loop
	    B (I) := B (I) + A (I + Offset_1, J + Offset_2) * X (J);
	 end loop;
      end loop;
      return B;
   end "*";
   
   function "*" (A : in Real_Matrix;
		 B : in Real_Matrix) return Real_Matrix is
      C : Real_Matrix (1 .. A'Length (1), 1 .. B'Length (2))
	:= (others => (others => 0.0));
      Aoff_1 : constant Integer := A'First (1) - 1;
      Boff_1 : constant Integer := B'First (1) - A'First (2);
      Boff_2 : constant Integer := B'First (2) - 1;
   begin
      for I in C'Range (1) loop
	 for J in C'Range (2) loop
	    for K in A'Range (2) loop
	       C (I, J) := C (I, J) 
		 + A (I + Aoff_1, K) * B (K + Boff_1, J + Boff_2);
	    end loop;
	 end loop;
      end loop;
      return C;
   end "*";


   function Transpose (A : in Real_Matrix) return Real_Matrix is
      B : Real_Matrix (A'Range (2), A'Range (1));
   begin
      for I in A'Range (1) loop
	 for J in A'Range (2) loop
	    B (J, I) := A (I, J);
	 end loop;
      end loop;
      return B;
   end Transpose;
   
   function Eye (N : in Pos) return Real_Matrix is
      A : Real_Matrix (1 .. N, 1 .. N) := (others => (others => 0.0));
   begin
      for I in 1 .. N loop
	 A (I, I) := 1.0;
      end loop;
      return A;
   end Eye;
   
   function Eye (N : in Pos) return Int_Matrix is
      A : Int_Matrix (1 .. N, 1 .. N) := (others => (others => 0));
   begin
      for I in 1 .. N loop
	 A (I, I) := 1;
      end loop;
      return A;
   end Eye;
      
   
   function Pivoting_Array (A : in Real_Matrix) return Int_Array is
      N   : constant Nat := A'Length (1);
      P   : Int_Array (1 .. N);
      Max : Real;
      Row : Pos;
      Tmp : Nat;
   begin
      for I in P'Range loop
	 P (I) := I;
      end loop;
      
      for J in 1 .. N loop
         Max := A (J + A'First (1) - 1, J + A'First (2) - 1);
         Row := J;
         for I in J + 1 .. N loop
            if A (I + A'First (1) - 1, J + A'First (2) - 1) > Max then
               Max := A (I + A'First (1) - 1, J + A'First (2) - 1);
               Row := I;
            end if;
         end loop;
	 if J /= Row then
	    Tmp     := P (J);
	    P (J)   := P (Row);
	    P (Row) := Tmp;
         end if;
      end loop;
      return P;
   end Pivoting_Array;
   
   function Permute_Row (A : in Real_Matrix;
			 P : in Int_Array) return Real_Matrix is
      PA : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
   begin
      for Row in P'Range loop
	 for J in 1 .. A'Length (1) loop
	    PA (Row, J) := A (P (Row), J + A'First (1) - 1);
	 end loop;
      end loop;
      return PA;
   end Permute_Row;
   
   procedure LU_Decomposition (A : in     Real_Matrix;
			       P :    out Int_Array;
			       L :    out Real_Matrix;
			       U :    out Real_Matrix) is
      N  : constant Nat := A'Length (1);
      PA : Real_Matrix (A'Range (1), A'Range (2));
      S  : Real;
   begin
      L := (others => (others => 0.0));
      U := (others => (others => 0.0));
      P := Pivoting_Array (A);
      PA := Permute_Row (A => A, P => P);
      for J in 0 .. N - 1 loop
         L (L'First (1) + J, L'First (2) + J) := 1.0;
         for I in 0 .. J loop
            S := 0.0;
            for K in 0 .. I - 1 loop
               S := S + U (U'First (1) + K, U'First (2) + J) *
                 L (L'First (1) + I, L'First (2) + K);
            end loop;
            U (U'First (1) + I, U'First (2) + J) :=
              PA (PA'First (1) + I, PA'First (2) + J) - S;
         end loop;
         for I in J + 1 .. N - 1 loop
            S := 0.0;
            for K in 0 .. J loop
               S := S + U (U'First (1) + K, U'First (2) + J) *
                 L (L'First (1) + I, L'First (2) + K);
            end loop;
            L (L'First (1) + I, L'First (2) + J) :=
              (PA (PA'First (1) + I, PA'First (2) + J) - S) /
              U (U'First (1) + J, U'First (2) + J);
         end loop;
      end loop;
   end LU_Decomposition;
   
   procedure Print (A : in Real_Matrix) is
      use Real_IO, Ada.Text_IO;
   begin
      for I in A'Range (1) loop
	 for J in A'Range (2) loop
	    Put (A (I, J), Aft => 3, Exp => 0); Put (",  ");
	 end loop;
	 New_Line;
      end loop;
   end Print;
   
   procedure Print (A : in Int_Matrix) is
      use Int_IO, Ada.Text_IO;
   begin
      for I in A'Range (1) loop
	 for J in A'Range (2) loop
	    Put (A (I, J), 3); Put (",  ");
	 end loop;
	 New_Line;
      end loop;
   end Print;
   
   function Number_Of_Swaps (P : in Int_Array) return Pos is
      Y   : array (P'Range) of Boolean := (others => False);
      Ind : Pos := 0;
      Cycles : Pos := 0;
   begin
      for I in P'Range loop
	 if Y (I) = False then
	    Cycles := Cycles + 1;
	    Y (I)  := True;
	    Ind    := P (I);
	    while Y (Ind) = False loop
	       Y (Ind) := True;
	       Ind := P (Ind);
	    end loop;
	 end if;
      end loop;
      return P'Length - Cycles;
   end Number_Of_Swaps;
   
   function Diag (A : in Real_Matrix) return Real_Vector is
      X : Real_Vector (1 .. A'Length (1));
   begin
      for I in X'Range loop
	 X (I) := A (I + A'First (1) - 1, I + A'First (2) - 1);
      end loop;
      return X;
   end Diag;
   
   function Diag (X : in Real_Vector) return Real_Matrix is
      A : Real_Matrix (1 .. X'Length, 1 .. X'Length)
	:= (others => (others => 0.0));
   begin
      for I in X'Range loop
	 A (I + 1 - X'First, I + 1 - X'First) := X (I);
      end loop;
      return A;
   end Diag;
   
   function Determinant (P : in Int_Array;
			 L : in Real_Matrix;
			 U : in Real_Matrix) return Real is
      Num : Pos := Number_Of_Swaps (P);
      Vec_L : Real_Vector := Diag (L);
      Vec_U : Real_Vector := Diag (U);
      Det : Real := 1.0;
   begin
      for X of Vec_L loop
	 Det := Det * X;
      end loop;
      for X of Vec_U loop
	 Det := Det * X;
      end loop;
      if Num mod 2 = 1 then Det := -Det; end if;
      return Det;
   end Determinant;
   
   function Determinant (A : in Real_Matrix) return Real is
      P : Int_Array (1 .. A'Length (1));
      L : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
      U : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
   begin
      LU_Decomposition (A, P, L, U);
      return Determinant (P, L, U);
   end Determinant;
   
   function Solve_Upper (U : in Real_Matrix;
			 B : in Real_Vector) return Real_Vector is
      N : constant Nat := U'Length (1);
      X : Real_Vector (1 .. N) := (others => 0.0);
      Tmp : Real;
   begin
      
      for I in reverse 1 .. N loop
	 Tmp := 0.0;
	 for J in I + 1 .. N loop
	    Tmp := Tmp + U (I, J) * X (J);
	 end loop;
	 X (I) := (B (I) - Tmp) / U (I, I);
      end loop;
      return X;
   end Solve_Upper;
   
   function Solve_Lower (L : in Real_Matrix;
			 B : in Real_Vector) return Real_Vector is
      N : constant Nat := L'Length (1);
      X : Real_Vector (1 .. N) := (others => 0.0);
      Tmp : Real;
   begin
      for I in 1 .. N loop
	 Tmp := 0.0;
	 for J in 1 .. I - 1 loop
	    Tmp := Tmp + L (I, J) * X (J);
	 end loop;
	 X (I) := B (I) - Tmp;
      end loop;
      return X;
   end Solve_Lower;
   
   function Permute_Row (X : in Real_Vector;
			 P : in Int_Array) return Real_Vector is
      Y : Real_Vector (X'Range);
   begin
      for I in P'Range loop
	 Y (I) := X (P (I));
      end loop;
      return Y;
   end Permute_Row;
   
   function Solve (P : in Int_Array;
   		   L : in Real_Matrix;
   		   U : in Real_Matrix;
   		   B : in Real_Vector) return Real_Vector is
      Pb : Real_Vector := Permute_Row (B, P);
      Y  : Real_Vector (1 .. B'Length);
   begin      
      Y := Solve_Lower (L, Pb);
      return Solve_Upper (U, Y);
   end Solve;
   
   function Solve (A : in Real_Matrix;
		   B : in Real_Vector) return Real_Vector is
      P : Int_Array (1 .. A'Length (1));
      L : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
      U : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
   begin
      LU_Decomposition (A, P, L, U);
      return Solve (P, L, U, B);
   end Solve;
   
   function Inverse (A : in Real_Matrix) return Real_Matrix is
      P : Int_Array (1 .. A'Length (1));
      X : Real_Vector (1 .. A'Length (1));
      L : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
      U : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
      Inv : Real_Matrix (1 .. A'Length (1), 1 .. A'Length (2));
   begin
      LU_Decomposition (A, P, L, U);
      pragma Assert (Determinant (P, L, U) /= 0.0, "Error: matrix is singular");
      
      for I in X'Range loop
	 for Item of X loop Item := 0.0; end loop;
	 X (I) := 1.0;
	 X := Solve (P, L, U, X);
	 for J in X'Range loop
	    Inv (J, I) := X (J);
	 end loop;
      end loop;
      return Inv;
   end Inverse;
end Numerics.Dense_Matrices;
