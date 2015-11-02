with Ada.Numerics.Generic_Real_Arrays;
package body Numerics is
   
   package Real_Arrays is new Ada.Numerics.Generic_Real_Arrays (Real);
   
   function Rand return Real is 
      use Ada.Numerics.Float_Random;
   begin
      return Real (Random (Gen));
   end Rand;
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   procedure Set (X  : in out Real_Vector;
		  To : in     Real_Array) is
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   procedure Set (X  : in out Int_Vector;
		  To : in     Int_Array) is
   begin
      X.Reserve_Capacity (To'Length);
      X.Set_Length (0);
      for Y of To loop
	 X.Append (Y);
      end loop;
   end Set;
   
   function Vectorize (Item : in Real_Array) return Real_Vector is
      X : Real_Vector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return Int_Vector is
      X : Int_Vector;
   begin
      X.Reserve_Capacity (Item'Length);
      X.Set_Length (0);
      for Y of Item loop
	 X.Append (Y);
      end loop;
      return X;
   end Vectorize;
   
   function Norm (X : in Real_Array) return Real is
      Y : Real := 0.0;
      use Real_Functions;
   begin
      for Item of X loop
	 Y := Y + Item ** 2;
      end loop;
      return Sqrt (Y);
   end Norm;

   
   function To_Array (Item : in Real_Vector) return Real_Array is
      Result : Real_Array (1 .. Nat (Item.Length));
      I : Nat := 1;
      use Ada.Text_IO;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   function To_Array (Item : in Int_Vector) return Int_Array is
      Result : Int_Array (1 .. Nat (Item.Length));
      I : Nat := 1;
      use Ada.Text_IO;
   begin
      for Y of Item loop
	 Result (I) := Y;
	 I := I + 1;
      end loop;
      return Result;
   end To_Array;
   
   
   
   ----- Vector and Array functions
   
   
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int) is
   begin
      V.Set_Length (Ada.Containers.Count_Type (N));
   end Set_Length;

   function Length (X : in Real_Vector) return Int is (Int (X.Length));
   
   
   
   -------- Max and Abs_Max functions ------------------
   
   function Max_Int_Array (Item : in Int_Array) return Int is separate;
   function Max_Real_Array (Item : in Real_Array) return Real is separate;
   function Abs_Max_IA (Item : in Int_Array) return Int is separate;
   function Abs_Max_RA (Item : in Real_Array) return Real is separate;
   
   function Max (X : in Int_Vector) return Int is
      Result : Int;
   begin
      case X.Length is
   	 when 0 => return 0;
   	 when 1 => return X (1);
   	 when others =>
   	    Result := X (1);
   	    for I in 2 .. Int (X.Length) loop
   	       if X (I) > Result then Result := X (I); end if;
   	    end loop;
   	    return Result;
      end case;
   end Max;
   
   function Max (X : in Real_Vector) return Real is
      Result : Real;
   begin
      case X.Length is
   	 when 0 => return 0.0;
   	 when 1 => return X (1);
   	 when others =>
   	    Result := X (1);
   	    for I in 2 .. Int (X.Length) loop
   	       if X (I) > Result then Result := X (I); end if;
   	    end loop;
   	    return Result;
      end case;
   end Max;


   -------- Binary Operators ---------------------------
   function Dot_Product (Left_I, Right_J : in Int_Array;
   			 Left_X, Right_Y : in Real_Array) return Real is separate;
   
   
   function Sparse (X	: in Real_Array;
		    N	: in Pos	:= 0;
		    Tol	: in Real	:= 1.0e-20) return Sparse_Vector is
      use IV_Package, RV_Package, Ada.Containers;
      Y : Sparse_Vector;
   begin
      Y.NMax := (if N < X'Length then X'Length else N);
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
   
   
   procedure Print (X : in Sparse_Vector) is
      use Int_IO, Real_IO, Ada.Text_IO;
   begin
      Put ("Length of vector:"); Put (X.NMax); New_Line;
      for I in 1 .. Pos (X.X.Length) loop
	 Put (X.I (I)); Put (", "); Put (X.X (I)); New_Line;
      end loop;
   end Print;
   
   function Length (X : in Sparse_Vector) return Pos is (X.NMax);

   
   function To_Array (X	  : in Sparse_Vector) return Real_Array is
      Y : Real_Array (1 .. X.NMax);
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
      elsif Item.I.Contains (I) = False then
	 for K in 1 .. Length loop
	    if Item.I (K) > I then
	       Item.I.Insert (Before => K, New_Item => I);
	       Item.X.Insert (Before => K, New_Item => X);
	       return;
	    end if;
	 end loop;
	 Item.I.Append (I);
	 Item.X.Append (X);
      else
	 Item.X (Item.I.Find_Index (I)) := X;
      end if;
   end Set;
   

begin
   Ada.Numerics.Float_Random.Reset (Gen);
end Numerics;
