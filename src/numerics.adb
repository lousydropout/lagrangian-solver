package body Numerics is
   
   function Rand return Real is 
      use Ada.Numerics.Float_Random;
   begin
      return Real (Random (Gen));
   end Rand;
   
   function Rand (N : in Nat) return Real_Vector is
      use Ada.Numerics.Float_Random;
      X : Real_Vector (1 .. N) := (others => Rand);
   begin
      return X;
   end Rand;
   
   
   --  function Solve (A : in Real_Matrix;
   --  		   B : in Real_Vector) return Real_Vector is
   --     C : RA.Real_Vector := RA.Real_Vector (B);
   --  begin
   --     C := RA.Solve (A, C);
   --     return Real_Vector (C);
   --  end Solve;
   
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out RVector;
		  To : in     Real_Vector) is
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   procedure Set (X  : in out IVector;
		  To : in     Int_Array) is
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   function Vectorize (Item : in Real_Vector) return RVector is
      X : RVector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return IVector is
      X : IVector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function Norm (X : in Real_Vector) return Real is
      Y : Real := 0.0;
      use Real_Functions;
   begin
      for Item of X loop
	 Y := Y + Item ** 2;
      end loop;
      return Sqrt (Y);
   end Norm;

   
   function To_Array (Item : in RVector) return Real_Vector is
      Result : Real_Vector (1 .. Pos (Item.Length));
      I : Nat := 1;
      use Ada.Text_IO;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   function To_Array (Item : in IVector) return Int_Array is
      Result : Int_Array (1 .. Pos (Item.Length));
      I : Nat := 1;
      use Ada.Text_IO;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   function Zero (N : in Pos) return Sparse_Vector is
      X : Real_Vector (1 .. 1) := (others => 0.0);
      Y : Sparse_Vector := Sparse (X);
   begin
      Y.NMax := N;
      return Y;
   end Zero;
   
   function One (I, N : in Nat) return Sparse_Vector is
      X : Real_Vector (1 .. 1) := (others => 1.0);
      Y : Sparse_Vector;
   begin
      Y := Zero (I - 1) or Sparse (X) or Zero (N - I);
      return Y;
   end One;
   ----- Vector and Array functions
   
   
   procedure Set_Length (V : in out RVector;
			 N : in     Integer) is
   begin
      V.Set_Length (Ada.Containers.Count_Type (N));
   end Set_Length;

   function Length (X : in RVector) return Integer is (Integer (X.Length));
   
   
   
   -------- Max and Abs_Max functions ------------------
   
   function Max_Int_Array (Item : in Int_Array) return Integer is separate;
   function Max_Real_Array (Item : in Real_Vector) return Real is separate;
   function Abs_Max_IA (Item : in Int_Array) return Integer is separate;
   function Abs_Max_RA (Item : in Real_Vector) return Real is separate;
   
   function Max (X : in IVector) return Integer is
      Result : Integer;
   begin
      case X.Length is
   	 when 0 => return 0;
   	 when 1 => return X (1);
   	 when others =>
   	    Result := X (1);
   	    for I in 2 .. Nat (X.Length) loop
   	       if X (I) > Result then Result := X (I); end if;
   	    end loop;
   	    return Result;
      end case;
   end Max;
   
   function Max (X : in RVector) return Real is
      Result : Real;
   begin
      case X.Length is
   	 when 0 => return 0.0;
   	 when 1 => return X (1);
   	 when others =>
   	    Result := X (1);
   	    for I in 2 .. Nat (X.Length) loop
   	       if X (I) > Result then Result := X (I); end if;
   	    end loop;
   	    return Result;
      end case;
   end Max;


   -------- Binary Operators ---------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
   			 Left_X, Right_Y : in Real_Vector) return Real is separate;
   
   
   function Sparse (X	: in Real_Vector;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 10.0 * Real'Small)
		   return Sparse_Vector is
      use IV_Package, RV_Package, Ada.Containers;
      Y : Sparse_Vector;
      Offset : constant Integer := 1 - X'First;
   begin
      Y.NMax := (if N < X'Length then X'Length else N);
      Y.X.Reserve_Capacity (Count_Type (X'Length));
      Y.I.Reserve_Capacity (Count_Type (X'Length));
      for I in X'Range loop
	 if abs (X (I)) > Tol then
	    Y.X.Append (X (I));
	    Y.I.Append (I + Offset);
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
      N  : constant Count_Type := Count_Type'Min (Count_Type (A.NMax),
						  A.X.Length + B.X.Length);
   begin
      pragma Assert (A.NMax = B.NMax,
		     "ERROR: Vectors are not of equal lengths");
      C.NMax := A.NMax;
      C.X.Reserve_Capacity (N);
      C.I.Reserve_Capacity (N);
      
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
   
   
   function "*" (A : in Real_Matrix;
		 B : in Sparse_Vector) return Sparse_Vector is
      C : Real_Vector (A'Range (2)) := (others => 0.0);
      J : constant Int_Array  := To_Array (B.I);
      X : constant Real_Vector := To_Array (B.X);
   begin
      
      for K in J'Range loop
	 for I in A'Range (1) loop
	    C (I) := C (I) + A (I, J (K)) * X (J (K));
	 end loop;
      end loop;
      
      return Sparse (C);
   end "*";

   
   
   procedure Print (X : in Sparse_Vector) is
      use Int_IO, Real_IO, Ada.Text_IO;
   begin
      Put ("Length of vector:"); Put (X.NMax); New_Line;
      for I in 1 .. Pos (X.X.Length) loop
	 Put (X.I (I)); Put (", "); Put (X.X (I)); New_Line;
      end loop;
   end Print;
   
   function Length (X : in Sparse_Vector) return Pos is (X.NMax);

   
   function To_Array (X	  : in Sparse_Vector) return Real_Vector is
      Y : Real_Vector (1 .. X.NMax);
   begin
      if Pos (X.X.Length) = X.NMax then
	 Y := To_Array (X.X);
      else
	 Y := (others => 0.0);
	 for K in 1 .. Pos (X.I.Length) loop
	    Y (X.I (K)) := X.X (K);
	 end loop;
      end if;

      return Y;
   end To_Array;
   
   function Norm (X : in Sparse_Vector) return Real is
      Result : Real := 0.0;
   begin
      for Item of X.X loop
	 Result := Result + Item ** 2;
      end loop;
      return Real_Functions.Sqrt (Result);
   end Norm;
   
   procedure Set_Length (X : in out Sparse_Vector;
			 N : in     Pos) is
      use Ada.Containers;
   begin
      X.NMax := N;
      X.X.Reserve_Capacity (Count_Type (N));
      X.I.Reserve_Capacity (Count_Type (N));
   end Set_Length;
   
   
   procedure Set (Item : in out Sparse_Vector;
   		  I    : in     Nat;
   		  X    : in     Real) is
      Length : constant Pos := Pos (Item.I.Length);
   begin
      if Length = Item.NMax then
	 Item.X (I) := X;
      elsif Length = 0 then
	 Item.I.Append (I);
	 Item.X.Append (X);
      else
	 -- look for index K >= I in Item.I
	 for K in 1 .. Length loop
	    if Item.I (K) = I then
	       Item.X (K) := X;
	       return;
	    elsif Item.I (K) > I then
	       Item.I.Insert (Before => K, New_Item => I);
	       Item.X.Insert (Before => K, New_Item => X);
	       return;
	    end if;
	 end loop;
	 -- if not found, then append since I is largest
	 Item.I.Append (I);
	 Item.X.Append (X);
      end if;
   end Set;
   

   procedure Add (Item : in out Sparse_Vector;
   		  I    : in     Nat;
   		  X    : in     Real) is
      Length : constant Pos := Pos (Item.I.Length);
   begin
      if Length = Item.NMax then
	 Item.X (I) := X;
      elsif Length = 0 then
	 Item.I.Append (I);
	 Item.X.Append (X);
      else
	 -- look for index K >= I in Item.I
	 for K in 1 .. Length loop
	    if Item.I (K) = I then
	       Item.X (K) := Item.X (K) + X;
	       return;
	    elsif Item.I (K) > I then
	       Item.I.Insert (Before => K, New_Item => I);
	       Item.X.Insert (Before => K, New_Item => X);
	       return;
	    end if;
	 end loop;
	 Item.I.Append (I);
	 Item.X.Append (X);
      end if;
   end Add;
   
   
   function Remove_1st_N (X : in Sparse_Vector;
			  N : in Pos) return Sparse_Vector is
      Y : Sparse_Vector;
   begin
      pragma Assert (X.NMax > N);
      Y.NMax := X.NMax - N;
      Y.I.Reserve_Capacity (X.X.Length);
      Y.X.Reserve_Capacity (X.X.Length);
      
      for I in 1 .. Pos (X.X.Length) loop
	 if X.I (I) > N then
	    Y.I.Append (X.I (I) - N);
	    Y.X.Append (X.X (I));
	 end if;
      end loop;
      Y.I.Reserve_Capacity (Y.I.Length);
      Y.X.Reserve_Capacity (Y.X.Length);
      return Y;
   end Remove_1st_N;
   
   function "or" (X, Y : in Sparse_Vector) return Sparse_Vector is
      use Ada.Containers;
      Z : Sparse_Vector;
   begin
      Z.NMax := X.NMax + Y.NMax;
      Z.I.Reserve_Capacity (X.I.Length + Y.I.Length);
      Z.X.Reserve_Capacity (X.I.Length + Y.I.Length);
      
      for I in 1 .. Pos (X.X.Length) loop
	 Z.I.Append (X.I (I));
	 Z.X.Append (X.X (I));
      end loop;
      for I in 1 .. Pos (Y.X.Length) loop
	 Z.I.Append (Y.I (I) + X.NMax);
	 Z.X.Append (Y.X (I));
      end loop;
      return Z;
   end "or";
   
   function "+" (X, Y : in Pos2D) return Pos2D is
   begin
      return (X.X + Y.X, X.Y + Y.Y);
   end "+";
   
   function "-" (X : in Pos2D) return Pos2D is
   begin
      return (-X.X, -X.Y);
   end "-";
   
   function "-" (X, Y : in Pos2D) return Pos2D is
   begin
      return (X.X - Y.X, X.Y - Y.Y);
   end "-";

   function Dot (X, Y : in Pos2D) return Real is
   begin
      return X.X * Y.X + X.Y * Y.Y;
   end Dot;
   
   function Norm (X : in Pos2D) return Real is
      use Real_Functions;
   begin
      return Sqrt (X.X**2 + X.Y**2);
   end Norm;
   
   function To_Array (Xvec : in Pos2D_Vector) return Real_Vector is
      Result : Real_Vector (1 .. 2 * Xvec'Length);
      K      : Nat := 1;
   begin
      for X of Xvec loop
	 Result (K) := X.X; K := K + 1;
	 Result (K) := X.Y; K := K + 1;
      end loop;
      return Result;
   end To_Array;
   
   
   function "+" (X, Y : in Pos2D_Vector) return Pos2D_Vector is
      Z : Pos2D_Vector := X;
   begin
      for I in Z'Range loop
	 Z (I) := Z (I) + Y (I);
      end loop;
      return Z;
   end "+";
   
   function "-" (X, Y : in Pos2D_Vector) return Pos2D_Vector is
      Z : Pos2D_Vector := X;
   begin
      for I in Z'Range loop
	 Z (I) := Z (I) - Y (I);
      end loop;
      return Z;
   end "-";
   
   function "-" (X : in Pos2D_Vector) return Pos2D_Vector is
      Z : Pos2D_Vector (X'Range);
   begin
      for I in Z'Range loop
	 Z (I) := -X (I);
      end loop;
      return Z;
   end "-";
   
   
   
   
   function "*" (X : in Real;
		 Y : in Pos2D_Vector) return Pos2D_Vector is
      Z : Pos2D_Vector := Y;
   begin
      for K of Z loop
	 K := X * K;
      end loop;
      return Z;
   end "*";
   
   
   -------- Real array functions ----------------
   function "-" (X : in Real_Vector) return Real_Vector is
      Y : Real_Vector (X'Range);
   begin
      for I in X'Range loop
	 Y (I) := -X (I);
      end loop;
      return Y;
   end "-";
   
   function "+" (X : in Real_Vector;
		 Y : in Real_Vector) return Real_Vector is
      Z : Real_Vector (1 .. X'Length);
      Offset_Y : constant Integer := Y'First - X'First;
      Offset_Z : constant Integer := 1 - X'First;
   begin
      for I in X'Range loop
	 Z (I + Offset_Z) := X (I) + Y (I + Offset_Y);
      end loop;
      return Z;
   end "+";
   
   function "-" (X : in Real_Vector;
		 Y : in Real_Vector) return Real_Vector is
      Z : Real_Vector (1 .. X'Length);
      Offset_Y : constant Integer := Y'First - X'First;
      Offset_Z : constant Integer := 1 - X'First;
   begin
      for I in X'Range loop
	 Z (I + Offset_Z) := X (I) - Y (I + Offset_Y);
      end loop;
      return Z;
   end "-";
   
   function "*" (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector is
      Result : Real_Vector (Right'Range);
   begin
      for I in Result'Range loop
	 Result (I) := Left * Right (I);
      end loop;
      return Result;
   end "*";

   -----------   Real_Vector functions ---------------------------
   function "*" (A : in Real_Matrix;
		 X : in Real_Vector) return Real_Vector is
      Y : Real_Vector (1 .. A'Length (1)) := (others => 0.0);
      Offset1 : constant Integer := A'First (1) - 1;
      Offset2 : constant Integer := A'First (2) - X'First;
   begin
      for I in Y'Range loop
	 for J in X'Range loop
	    Y (I) := Y (I) + A (I + Offset1, J + Offset2) * X (J);
	 end loop;
      end loop;
      return Y;
   end "*";
   
   function Dot (X, Y : in Real_Vector) return Real is
      R : Real := 0.0;
   begin
      for I in X'Range loop
	 R := R + X (I) * Y (I + Y'First - X'First);
      end loop;
      return R;
   end Dot;
   
   
   
begin
   Ada.Numerics.Float_Random.Reset (Gen);
end Numerics;
