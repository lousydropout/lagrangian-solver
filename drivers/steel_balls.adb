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
   A, Q : Array_Of_Vectors;
   File : File_Type;
   Fcsv : File_Type;
   Fp   : File_Type;
   Dt   : Real;
   Name : String := "out.xyz";
   Cname : String := "out.csv";
   Total_Energy, T_Final  : Real;
   Line : String (1 .. 50);
   Last : Natural;
   Upper, Lower : Real;
   
begin
   Control.Max_Dt := 1.0e2;
   Control.Tol    := 1.0e-12;
   ------------------------------------------------------------
   -- Read initial conditions
   T     := 0.0;
   X (1) := 0.0;
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- s
   X (2) := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- ω_s
   X (4) := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- total energy
   Total_Energy := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- alpha
   α  := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- dt
   Dt := Real'Value (Line (1 .. Last));
   
   Get_Line (Line, Last); -- Not used
   Get_Line (Line, Last); -- T_final
   T_Final := Real'Value (Line (1 .. Last));
   ------------------------------------------------------------
   X := Get_IC (X, Total_Energy);
   Total_Energy := Val (Hamiltonian (0.0, X));
   Put ("Total Energy = "); Put (Total_Energy); New_Line;
   ------------------------------------------------------------
   State   := Var;
   ------------------------------------------------------------
   ------------------------------------------------------------
   Create (Fp, Name => "poincare.csv");
   Put_Line (Fp, "time, t, s, t_dot, s_dot, pt, ps, E");
   Print_Lagrangian (Fp, Var, Lagrangian'Access);
   Print_XYZ (File, State, Name, Create);
   Print_CSV (Fcsv, State, CName, Lagrangian'Access, Create);
   
   while T < T_Final loop
      Y := Update (Lagrangian'Access, Var, Control, Sparse);
      A := Chebyshev_Transform (Y);
      Put (Var.T); New_Line;
      ---------------------------------------------------------
      Q := Split (Y);
      for I in 2 .. K loop
   	 if X1 (Q (1) (I - 1)) * X1 (Q (1) (I)) < 0.0 and then 
   	   Y1 (Q (1) (I)) > 0.0 then
   	 -- If there's a zero, bisect
   	    Put ("****************************    ");
   	    Lower   := Var.T + Control.Dt * Grid (I - 1);
   	    Upper   := Var.T + Control.Dt * Grid (I);
   	    Guess   := Find_State_At_Level
   	      (0.0, A, Var.T, Control.Dt, Lower, Upper, X1'Access);
   	    ----------------------------------------------------------------
   	    Put (Guess.X (3)); New_Line;
   	    if Guess.X (3) > 0.0 then
   	       Print_Lagrangian (Fp, Guess, Lagrangian'Access);
   	    end if;
      	 end if;
      end loop;
      ---------------------------------------------------------
      
      
      ---------------------------------------------------------
      while State.T + Dt <= T + Control.Dt loop
      	 State.T := State.T + Dt;
   	 State.X := Interpolate (A, State.T, Var.T, Var.T + Control.Dt);
   	 Print_XYZ (File, State, Name);
   	 Print_CSV (Fcsv, State, CName, Lagrangian'Access);
      end loop;
      ---------------------------------------------------------
      
      Update (Var => Var, Y => Y, Dt => Control.Dt); -- Update variable Var
   end loop;
    
   Close (Fp);
end Steel_Balls;
