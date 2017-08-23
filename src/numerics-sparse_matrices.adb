package body Numerics.Sparse_Matrices is
   
   procedure Print (Mat : in Sparse_Matrix) is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Sparse_Matrix) return Real is separate;
   function N_Row (Mat : in Sparse_Matrix)  return Pos  is separate;
   function N_Col (Mat : in Sparse_Matrix)  return Pos  is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Functions for Creating Sparse Matrices -------------------
   function Sparse (X : in Real_Matrix) return Sparse_Matrix is
      use Ada.Containers;
      Y : Sparse_Matrix;
      N : constant Count_Type := Count_Type (X'Length (1) * X'Length (2));
   begin
      Y.N_Row := X'Length (1);
      Y.N_Col := X'Length (2);
      Y.Format := Triplet;
      Y.I.Reserve_Capacity (N);
      Y.P.Reserve_Capacity (N);
      Y.X.Reserve_Capacity (N);
      for I in 1 .. Int (X'Length (1)) loop
	 for J in 1 .. Int (X'Length (2)) loop
	    Y.I.Append (I); 
	    Y.P.Append (J);
	    Y.X.Append (X (I, J));
	 end loop;
      end loop;
      Y.Compress;
      return Y;
   end Sparse;
   function As_Matrix (X : in Sparse_Vector) return Sparse_Matrix is
      A : Sparse_Matrix;
      Eps : constant Real := 10.0 * Real'Small;
   begin
      A.Format := CSC; A.N_Col := 1; A.N_Row := X.NMax;
      A.X := X.X;
      A.I := X.I;
      A.P.Reserve_Capacity (2);
      A.P.Append (1);
      A.P.Append (Nat (X.I.Length) + 1);
      return A;
   end As_Matrix;
   
   function "*" (Left, Right : in Sparse_Vector) return Sparse_Matrix is
      use Ada.Containers;
      A, B : Sparse_Matrix;
      Eps : constant Real := 10.0 * Real'Small;
   begin
      A.Format := CSC; A.N_Col := 1; A.N_Row := Left.NMax;
      A.X := Left.X;
      A.I := Left.I;
      A.P.Reserve_Capacity (2);
      A.P.Append (1);
      A.P.Append (Nat (Left.I.Length) + 1);
	 
      
      B.Format := CSC; B.N_Col := 1; B.N_Row := Right.NMax;
      B.X := Right.X;
      B.I := Right.I;
      B.P.Reserve_Capacity (2);
      B.P.Append (1);
      B.P.Append (Nat (Right.I.Length) + 1);
      
      return (A * Transpose (B));
   end "*";
   

   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Array;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Sparse_Matrix_Format := CSC) 
			      return Sparse_Matrix is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   procedure Remove_Duplicates (Mat : in out Sparse_Matrix) is separate;
   procedure Compress (Mat : in out Sparse_Matrix) is separate;
   procedure Convert (Mat : in out Sparse_Matrix) is separate;
   function Convert (Mat : in Sparse_Matrix) return Sparse_Matrix is
      Result : Sparse_Matrix := Mat;
   begin
      Result.Convert;
      return Result;
   end Convert;
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Testing Functions -----------------------------------
   function Is_Square_Matrix (A : in Sparse_Matrix) return Boolean is separate;
   function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean is separate;   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix Operations -----------------------------------
   function Eye (N : in Nat) return Sparse_Matrix is separate;
   function Zero (N : in Nat) return Sparse_Matrix is separate;
   function Omega (N : in Nat;
		   M : in Pos := 0) return Sparse_Matrix is separate;

   procedure Transposed (Mat : in out Sparse_Matrix) is separate;
   function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Plus (Left  : in Sparse_Matrix;
		  Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Minus (Left  : in Sparse_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Direct_Sum (A, B : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Permute_By_Col (Mat : in Sparse_Matrix;
			    P   : in Int_Array) return Sparse_Matrix is separate;
   function Permute (Mat : in Sparse_Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Sparse_Matrix is separate;
   
   
   --  function BiCGSTAB (A   : in     Sparse_Matrix;
   --  		      B   : in     Real_Vector;
   --  		      X0  : in     Real_Vector;
   --  		      Err :    out Real;
   --  		      Tol : in     Real	    := 1.0e-10) 
   --  		     return Real_Vector is separate;

   function Number_Of_Elements (X : in Sparse_Matrix) return Int is (Int (X.X.Length));
   
   function Is_Valid (Mat : in Sparse_Matrix) return Boolean is
      use IV_Package, RV_Package;
   begin
      if Mat.I = IV_Package.Empty_Vector 
	or else Mat.P = IV_Package.Empty_Vector
	or else Mat.X = RV_Package.Empty_Vector 
	or else Mat.N_Row = 0 or else Mat.N_Col = 0 
      then
	 return False;
      end if;
      return True;
   end Is_Valid;
   
   procedure Triplet_To_Matrix (Result :    out Sparse_Matrix;
				I      : in     Int_Vector;
				J      : in     Int_Vector;
				X      : in     Real_Vector;
				N_Row  : in     Pos	      := 0;
				N_Col  : in     Pos	      := 0) is
   begin
      Result.N_Row  := Pos'Max (N_Row, Max (I));
      Result.N_Col  := Pos'Max (N_Col, Max (J));
      
      Result.Format := Triplet;
      Result.X := X; Result.I := I; Result.P := J; 
      
      Result.Compress; -- turns it into CSC format
   end Triplet_To_Matrix;
   
   function Read_Sparse_Triplet (File_Name : in String;
				 Offset	   : in Int    := 0)
				return Sparse_Matrix is
      use Ada.Text_IO, Ada.Containers, Real_IO, Int_IO;
      N_Lines : Count_Type := 0;
      I_Vec : Int_Vector;
      J_Vec : Int_Vector;
      X_Vec : Real_Vector;
      Int_Input : Int;
      Real_Input : Real;
      File : File_Type;
      Result : Sparse_Matrix;
   begin
      Open (File => File, Mode => In_File, Name => File_Name);

      while not End_Of_File (File) loop 
	 Get (File, Int_Input);  I_Vec.Append (Int_Input + 1 - Offset);
	 Get (File, Int_Input);  J_Vec.Append (Int_Input + 1 - Offset);
	 Get (File, Real_Input); X_Vec.Append (Real_Input);
	 N_Lines := N_Lines + 1;
      end loop;
      Close (File);
      
      I_Vec.Reserve_Capacity (N_Lines);
      J_Vec.Reserve_Capacity (N_Lines);
      X_Vec.Reserve_Capacity (N_Lines);
      
      Triplet_To_Matrix (Result, I_Vec, J_Vec, X_Vec);
      return Result;
   end Read_Sparse_Triplet;
   
   
   procedure Cumulative_Sum (Item : in out Int_Array) is separate;
   
   
   
   
   procedure Add (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real) is
      use Ada.Containers;
      Ind  : Pos;
   begin
      pragma Assert (Mat.Format = CSC);
      
      -- Check if Mat (I, J) exists
      for K in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (K) = I then 
	    -- If exists, then add X to Mat (I, J)
	    Mat.X (K) := Mat.X (K) + X;
	    return;
	 end if;
      end loop;
      
      -- Reserve space for 1 more element
      Mat.X.Reserve_Capacity (Mat.X.Length + 1);
      Mat.I.Reserve_Capacity (Mat.I.Length + 1);
      
      Ind := Mat.P (J); -- need this since for-loop may not occur
      for P in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (P) > I then Ind := P; exit; end if;
      end loop;
      
      -- Insert elements into I and X
      if Ind <= Pos (Mat.X.Length) then
	 Mat.X.Insert (Before => Ind, New_Item => X);
	 Mat.I.Insert (Before => Ind, New_Item => I);
      else
	 Mat.X.Append (X); Mat.I.Append (I);
      end if;

      -- Fix P
      for P in J + 1 .. Mat.N_Col + 1 loop
	 Mat.P (P) := Mat.P (P) + 1;
      end loop;      
      
   end Add;
   
   
   
   procedure Set (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real) is
      use Ada.Containers;
      Ind  : Pos;
   begin
      pragma Assert (Mat.Format = CSC);
      
      -- Check if Mat (I, J) exists
      for K in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (K) = I then 
	    -- If exists, then set Mat (I, J) to X
	    Mat.X (K) := X;
	    return;
	 end if;
      end loop;
      
      -- Reserve space for 1 more element
      Mat.X.Reserve_Capacity (Mat.X.Length + 1);
      Mat.I.Reserve_Capacity (Mat.I.Length + 1);
      
      Ind := Mat.P (J); -- need this since for-loop may not occur
      for P in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (P) > I then Ind := P; exit; end if;
      end loop;
      
      -- Insert elements into I and X
      if Ind <= Pos (Mat.X.Length) then
	 Mat.X.Insert (Before => Ind, New_Item => X);
	 Mat.I.Insert (Before => Ind, New_Item => I);
      else
	 Mat.X.Append (X); Mat.I.Append (I);
      end if;

      -- Fix P
      for P in J + 1 .. Mat.N_Col + 1 loop
	 Mat.P (P) := Mat.P (P) + 1;
      end loop;      
      
   end Set;
   
   
   
   procedure Scatter (A	   : in     Sparse_Matrix;
		      J	   : in     Int;
		      β	   : in     Real;
		      W	   : in out Int_Array;
		      X	   : in out Real_Array;
		      Mark : in     Int;
		      C	   : in out Sparse_Matrix;
		      Nz   : in out Int) is
      use IV_Package;
      I    : Int;
      Cur  : Cursor;
      L, R : Pos;
   begin
      Cur := To_Cursor (A.P, J);
      L   := A.P (Cur); Next (Cur); R := A.P (Cur) - 1;
      for P in L .. R loop
	 I := A.I (P);
	 if W (I) < Mark then
	    C.I.Append (I);
	    X (I) := β * A.X (P);
	    Nz    := Nz + 1;
	    W (I) := Mark;
	 else
	    X (I) := X (I) + β * A.X (P);
	 end if;
      end loop;
   end Scatter;
   
      
   function Mult_M_SV (A : in Sparse_Matrix;
		       X : in Sparse_Vector) return Sparse_Vector is separate;
   
   
   function Diag (X : in Sparse_Matrix) return Sparse_Vector is
      Y : Sparse_Vector;
   begin
      Set_Length (Y, X.N_Col);
      
      for P in 1 .. X.N_Col loop
	 for I in X.P (P) .. X.P (P + 1) - 1 loop
	    if X.I (I) = P then
	       Set (Y, P, X.X (I));
	       exit;
	    end if;
	 end loop;
      end loop;
      return Y;
   end Diag;
   
   function Diag (X : in Sparse_Vector) return Sparse_Matrix is
      Y : Sparse_Matrix;
   begin
      Y.Triplet_To_Matrix (X.I, X.I, X.X, X.NMax, X.NMax);
      return Y;
   end Diag;
      
   
   procedure Set_Diag (X  : in out Sparse_Matrix;
   		       To : in     Sparse_Vector) is
      K : Int;
   begin
      for I in 1 .. Int (To.I.Length) loop
   	 K := To.I (I);
   	 X.Set (K, K, To.X (I));
      end loop;
   end Set_Diag;
   
   
   function "-" (X : in Sparse_Matrix) return Sparse_Matrix is
      Y : Sparse_Matrix := X;
   begin
      for Item of Y.X loop
	 Item := -Item;
      end loop;
      return Y;
   end "-";
   
   function "*" (Left  : in Real;
		 Right : in Sparse_Matrix) return Sparse_Matrix is
      C : Sparse_Matrix := Right;
   begin
      for X of C.X loop
	 X := X * Left;
      end loop;
      return C;
   end "*";
   
   
   function Add_Column (X : in Sparse_Matrix;
			V : in Sparse_Vector) return Sparse_Matrix is
      use Ada.Containers;
      Y : Sparse_Matrix := X;
      N : Nat := Nat (X.X.Length);
      P : Nat := Nat (X.P.Length);
      M : Nat := Nat (V.X.Length);
   begin
      pragma Assert (X.N_Row = V.NMax);
      
      Y.N_Col := X.N_Col + 1;
      
      Y.X.Reserve_Capacity (Count_Type (N + V.NMax));
      Y.I.Reserve_Capacity (Count_Type (N + V.NMax));
      Y.P.Reserve_Capacity (Count_Type (P + 1));
      Y.P.Append (X.P (P) + M);
      
      for J in 1 .. M loop
	 Y.X.Append (V.X (J));
	 Y.I.Append (V.I (J));
      end loop;
      return Y;
   end Add_Column;
   
end Numerics.Sparse_Matrices;
