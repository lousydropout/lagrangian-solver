with Numerics, Ada.Text_IO, Chebyshev, Sb_Package;
use  Numerics, Ada.Text_IO, Chebyshev, Sb_Package;

procedure Steel_Balls is
   use Int_IO, Real_IO, AD_Package, Integrator;
   -----------------------------------------------
   Control : Control_Type := New_Control_Type;
   -- Initial Conditions ----
   Guess, Var, State : Variable;
   X : Vector renames Var.X;
   T : Real   renames Var.T;
   -------------------------------
   Y    : Real_Vector (1 .. NK);
   A    : Array_Of_Vectors;
   Fcsv : File_Type;
   Dt   : Real;
   Cname : String := "out.csv";
   Total_Energy, T_Final  : Real;
   Line : String (1 .. 50);
   Last : Natural;
   
begin
   Control.Max_Dt := 1.0e2;
   Control.Tol    := 1.0e-10;
   ------------------------------------------------------------
   -- Read initial conditions
   T     := 0.0;
   X (1) := 0.0;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- s
   X (2) := Real'Value (Line (1 .. Last));
   Put (X (2)); New_Line;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- ω_s
   X (4) := Real'Value (Line (1 .. Last));
   Put (X (4)); New_Line;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- total energy
   Total_Energy := Real'Value (Line (1 .. Last));
   Put ("Total_Energy = "); Put (Total_Energy); New_Line;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- alpha
   α  := Real'Value (Line (1 .. Last));
   Put ("alpha = "); Put (α); New_Line;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- dt
   Dt := Real'Value (Line (1 .. Last));
   Put (Dt); New_Line;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- T_final
   T_Final := Real'Value (Line (1 .. Last));
   Put ("T_Final = "); Put (T_Final); New_Line;
   
   ------------------------------------------------------------
   X := Get_IC (X, Total_Energy);
   --  X (3) := -X (3);
   for Item of X loop Put (Item); New_Line; end loop;
   Total_Energy := Val (Hamiltonian (0.0, X));
   Put ("Total Energy = "); Put (Total_Energy); New_Line;
   ------------------------------------------------------------
   State   := Var;
   ------------------------------------------------------------
   Print_CSV (Fcsv, State, CName, Lagrangian'Access, Create);
   
   while T < T_Final loop
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      A := Chebyshev_Transform (Y);
      Put (Var.T); Put ("             ");
      Put (Val (Hamiltonian (0.0, Var.X)));
      New_Line;
      ---------------------------------------------------------
      while State.T + Dt <= T + Control.Dt loop
      	 State.T := State.T + Dt;
   	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
   	 Print_CSV (Fcsv, State, CName, Lagrangian'Access);
      end loop;
      ---------------------------------------------------------
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
   
   Close (Fcsv);
end Steel_Balls;
