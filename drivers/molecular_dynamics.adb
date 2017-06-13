package body Molecular_Dynamics is
   
   function Verlet (R	  : in Pos2D_Vector;
		    R_New : in Pos2D_Vector;
		    Dt	  : in Real) return Pos2D_Vector is
      F : Pos2D_Vector := Calculate_Forces (R_New);
   begin
      return (2.0 * R_New - R + (Dt**2) * F);
   end Verlet;
   
   
   function LJ_Force (From : in Pos2D;
		      To   : in Pos2D) return Pos2D is
      use Real_Functions;
      E  : constant Real := 10000.0;
      R  : constant Real := Norm (To - From);
      R6 : constant Real := (1.0 / R) ** 6;
      F  : Real;
   begin
      F := (E / R ** 2) * R6 * (R6 - 1.0);
      return (F * (To - From));
   end LJ_Force;
   
   
   function Correct_Forces (F : in Pos2D_Vector;
			    M : in Nat) return Pos2D_Vector is
      G : Pos2D_Vector := F;
   begin
      G (1 .. 2 * M + 2) := (others => (0.0, 0.0));
      return G;
   end Correct_Forces;
   
   
   function Calculate_Forces (R : in Pos2D_Vector) return Pos2D_Vector is
      Force : Pos2D := (0.0, 0.0);
      F : Pos2D_Vector (R'Range) := (others => (0.0, 0.0));
      Check : Pos2D := (0.0, 0.0);
   begin
      for I in R'First .. R'Last loop
	 for J in I + 1 .. R'Last loop
	    Force := LJ_Force (From => R (I), To => R (J));
	    F (I) := F (I) - Force;
	    F (J) := F (J) + Force;
	 end loop;
      end loop;
      
      for X of F loop
	 Check := Check + X;
      end loop;
      if Norm (Check) > 0.1 then
	 Put_Line ("Error");
      end if;
      
      return Correct_Forces (F, 5);
   end Calculate_Forces;
   
   
   function Initialize_Lattice (N, M : in Nat;
				Vac  : in Nat) return Pos2D_Vector is
      K, Ind : Nat := 1;
      Result : Pos2D_Vector (1 .. 2 * N * M + N + M) := (others => (0.0, 0.0));
      Ratio : Real := 1.0;
   begin
      
      -- Row 2N
      for J in 0 .. M loop
	 Result (K).X := Ratio * Real (J);
	 Result (K).Y := Ratio * Real (2 * N) * Sqrt (3.0) * 0.5;
	 K := K + 1; Ind := Ind + 1;
      end loop;
      
      -- Row 0 
      for J in 0 .. M loop
	 Result (K).X := Ratio * Real (J);
	 Result (K).Y := 0.0;
	 K := K + 1; Ind := Ind + 1;
      end loop;
      
      -- The in-between rows
      for I in 1 .. 2 * N - 1 loop      
	 if I mod 2 = 0 then -- if row I is even
	    for J in 0 .. M loop
	       if Ind /= Vac then
		  Result (K).X := Ratio * Real (J);
		  Result (K).Y := Ratio * Real (I) * Sqrt (3.0) * 0.5;
		  K := K + 1; 
	       end if;
	       Ind := Ind + 1;
	    end loop;
	 else -- if row I is odd
	    for J in 1 .. M loop
	       if Ind /= Vac then
	    	  Result (K).X := Ratio * Real (J) - 0.5;
	    	  Result (K).Y := Ratio * Real (I) * Sqrt (3.0) * 0.5;
	    	  K := K + 1; 
	       end if;
	       Ind := Ind + 1;
	    end loop;
	    null;
	 end if;
      end loop;
      
      
      return Result;
   end Initialize_Lattice;
      

   
   
   
   
   
   
   
   procedure Output (R	  : in Pos2D_Vector;
		     File : in File_Type) is
   begin
      Put (File, R'Length); New_Line (File);
      Put_Line (File, "Properties=pos:R:2");
      
      for Item of R loop
	 Put (File => File, Item => Item.X);
	 Put (File => File, Item => "     ");
	 Put (File => File, Item => Item.Y);
	 Put (File => File, Item => "     ");
	 Put (File => File, Item => "0.5");
	 New_Line (File => File);
      end loop;
   end Output;
   

   
end Molecular_Dynamics;
