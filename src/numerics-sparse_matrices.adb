package body Numerics.Sparse_Matrices is
   
   procedure Print (Mat : in Sparse_Matrix) is separate;
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Sparse_Matrix) return Real is separate;
   function N_Row (Mat : in Sparse_Matrix) return Pos is separate;
   function N_Col (Mat : in Sparse_Matrix) return Pos is separate;
   
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
   --  function Cumulative_Sum (Item : in Int_Array) return Int_Array is separate;
   procedure Remove_Duplicates (Mat : in out Sparse_Matrix) is separate;
   procedure Compress (Mat : in out Sparse_Matrix) is separate;
   procedure Convert (Mat : in out Sparse_Matrix) is separate;
   function Convert (Mat : in Sparse_Matrix) return Sparse_Matrix is
      Result : Sparse_Matrix := Mat;
   begin
      Result.Convert;
      return Result;
   end Convert;
   
   
   
   function Vectorize (I : in Int_Array;
		       X : in Real_Array) return Sparse_Matrix is
      Result   : Sparse_Matrix;
      Offset_I : constant Int := I'First - 1;
      Offset_X : constant Int := X'First - 1;
   begin
      Result.Format := CSC;
      Result.N_Row := I (I'Last);
      Result.N_Col := 1;
      Result.P.Reserve_Capacity (2); 
      Result.I.Reserve_Capacity (I'Length);
      Result.X.Reserve_Capacity (X'Length);
      
      Result.P.Append (1);
      Result.P.Append (X'Length + 1);
      for Y of I loop
	 Result.I.Append (Y);
      end loop;
      for Y of X loop
	 Result.X.Append (Y);
      end loop;
      return Result;
   end Vectorize;

   
   
   
   
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   
   ------- Testing Functions -----------------------------------
   function Is_Col_Vector (A : in Sparse_Matrix) return Boolean is separate;
   function Is_Square_Matrix (A : in Sparse_Matrix) return Boolean is separate;
   function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean is separate;   
   
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix Operations -----------------------------------
   function Eye (N : in Nat) return Sparse_Matrix is separate;
   function Zero_Vector (N : in Nat) return Sparse_Matrix is separate;

   procedure Transposed (Mat : in out Sparse_Matrix) is separate;
   function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Plus (Left  : in Sparse_Matrix;
		  Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Minus (Left  : in Sparse_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Kronecker (A, B : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Direct_Sum (A, B : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Mult_M_RV (Left  : in Sparse_Matrix;
		       Right : in Real_Vector) return Real_Vector is separate;
   function Permute_By_Col (Mat : in Sparse_Matrix;
			    P   : in Int_Array) return Sparse_Matrix is separate;
   function Permute (Mat : in Sparse_Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Sparse_Matrix is separate;
   
   
   function BiCGSTAB (A   : in     Sparse_Matrix;
		      B   : in     Real_Vector;
		      X0  : in     Real_Vector;
		      Err :    out Real;
		      Tol : in     Real	    := 1.0e-10) 
		     return Real_Vector is separate;

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
   
   function Triplet_To_Matrix (I      : in Int_Vector;
			       J      : in Int_Vector;
			       X      : in Real_Vector;
			       N_Row  : in Pos		 := 0;
			       N_Col  : in Pos		 := 0;
			       Format : in Sparse_Matrix_Format := CSC) 
			      return Sparse_Matrix is
      Result : Sparse_Matrix;
   begin
      Result.N_Row  := (if N_Row = 0 then Max (I) else N_Row);
      Result.N_Col  := (if N_Col = 0 then Max (J) else N_Col);
      
      Result.Format := Triplet;
      Result.I := I; Result.P := J; Result.X := X;
      case Format is
	 when CSC     => 
	    Result.Compress;
	 when CSR     => 
	    Result.Compress; 
	    Result.Convert;
	 when Triplet => 
	    null;
      end case;
      return Result;
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
      
      return Triplet_To_Matrix (I_Vec, J_Vec, X_Vec);
   end Read_Sparse_Triplet;
   
   
   procedure Cumulative_Sum (Item : in out Int_Array) is separate;
   
   
   procedure Add (Mat  : in out Sparse_Matrix;
		  I, J : in     Nat;
		  X    : in     Real) is
      use RV_Package, IV_Package, Ada.Containers;
      Ind  : Pos;
      I1, I2 : Pos;
      X1, X2 : Real;
      CX : RV_Package.Cursor;
      CI : IV_Package.Cursor;
   begin
      pragma Assert (Mat.Format = CSC);
      -- Check if Mat (I, J) exists
      for K in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (K) = I then 
	    -- If exists, then just add value X to Mat (I, J)
	    Mat.X (K) := Mat.X (K) + X;
	    return;
	 end if;
      end loop;
      
      -- Fix P
      for P in J + 1 .. Mat.N_Col + 1 loop
	 Mat.P (P) := Mat.P (P) + 1;
      end loop;      
      
      -- Reserve space for 1 more element
      Mat.I.Reserve_Capacity (Mat.I.Length + 1);
      Mat.X.Reserve_Capacity (Mat.X.Length + 1);
      
      -- Fix I and X
      
      ---- 1. Find index of Mat.I & Mat.X just after (I, J)
      Ind := Mat.P (J);
      for P in Mat.P (J) .. Mat.P (J + 1) - 1 loop
	 if Mat.I (P) < I then Ind := P; end if;
      end loop;
      
      ---- 2. Get cursor for Mat.X at that index &
      --------- continuously swap value with previous value
      CX := To_Cursor (Mat.X, Ind);
      X2 := X;
      for K in Ind .. Nat (Mat.X.Length) loop
	 X1         := Mat.X (CX);
	 Mat.X (CX) := X2;
	 X2         := X1;
	 Next (CX);
      end loop;
      Mat.X.Append (X2);
      
      ---- 3. Repeat step 2 but for Mat.I
      CI := To_Cursor (Mat.I, Ind);
      I2 := I; 
      for K in Ind .. Nat (Mat.I.Length) loop
	 I1         := Mat.I (CI);
	 Mat.I (CI) := I2;
	 I2         := I1;
	 Next (CI);
      end loop;
      Mat.I.Append (I2);

   end Add;
   
   
   function Sparse (X	: in Real_Vector;
		    Tol	: in Real	 := 1.0e-10) return Sparse_Vector is
      use IV_Package, RV_Package, Ada.Text_IO;
      Y : Sparse_Vector;
   begin
      Put_Line ("Reserve cap");
      Y.X.Reserve_Capacity (X.Length);
      Y.I.Reserve_Capacity (X.Length);
      Put_Line ("begin loop");
      for I in 1 .. Nat (X.Length) loop
	 if abs (X (I)) > Tol then
	    Y.X.Append (X (I));
	    Y.I.Append (I);
	 end if;
      end loop;
      Put_Line ("finished");
      return Y;
   end Sparse;
   
   
   
   function Sparse (X	: in Real_Array;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-10) return Sparse_Vector is
      use IV_Package, RV_Package, Ada.Containers;
      Y : Sparse_Vector;
   begin
      Y.NMax := (if N = 0 then X'Length else N);
      Y.X.Reserve_Capacity (Count_Type (X'Length));
      Y.I.Reserve_Capacity (Count_Type (X'Length));
      for I in 1 .. Int (X'Length) loop
	 if abs (X (I)) > Tol then
	    Y.X.Append (X (I));
	    Y.I.Append (I);
	 end if;
      end loop;
      return Y;
   end Sparse;
   
   
   function "+" (A, B : in Sparse_Vector) return Sparse_Vector is
      use IV_Package, RV_Package, Ada.Containers;
      C : Sparse_Vector;
      Ax, Bx : Real;
      Ai, Bi : Pos;
      I, J : Pos := 1;
      Al : constant Pos := Pos (A.X.Length);
      Bl : constant Pos := Pos (B.X.Length);
   begin
      pragma Assert (A.NMax = B.NMax,
		     "ERROR: Vectors are not of equal lengths");
      C.NMax := A.NMax;
      C.X.Reserve_Capacity (A.X.Length + B.X.Length);
      C.I.Reserve_Capacity (A.X.Length + B.X.Length);
      
      while I <= Al and J <= Bl loop
	 Ax := A.X (I); Bx := B.X (J);
	 Ai := A.I (I); Bi := B.I (J);
	 
	 if Ai = Bi then
	    C.X.Append (Ax + Bx);
	    C.I.Append (Ai);
	    I := I + 1; J := J + 1;
	 elsif Bi < Ai then
	    C.X.Append (Bx);
	    C.I.Append (Bi);
	    J := J + 1;
	 else
	    C.X.Append (Ax);
	    C.I.Append (Ai);
	    I := I + 1;
	 end if;
      end loop;
      while I <= Al loop
	 C.X.Append (A.X (I));
	 C.I.Append (A.I (I));
	 I := I + 1;
      end loop;
      while J <= Bl loop
	 C.X.Append (B.X (J));
	 C.I.Append (B.I (J));
	 J := J + 1;
      end loop;
      C.X.Reserve_Capacity (C.X.Length);
      C.I.Reserve_Capacity (C.X.Length);
      
      return C;
   end "+";

   function "*" (A : in Real;
		 B : in Sparse_Vector) return Sparse_Vector is
      C : Sparse_Vector := B;
   begin
      for X of C.X loop
	 X := A * X;
      end loop;
      return C;
   end "*";
   
   
   procedure Print (X : in Sparse_Vector) is
      use Int_IO, Real_IO, Ada.Text_IO;
   begin
      Put ("Length of vector:"); Put (X.NMax); New_Line;
      for I in 1 .. Pos (X.X.Length) loop
	 Put (X.I (I)); Put (", "); Put (X.X (I)); New_Line;
      end loop;
   end Print;
end Numerics.Sparse_Matrices;
