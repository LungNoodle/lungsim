
.. |check| raw:: html

    <input checked=""  type="checkbox">

.. |check_| raw:: html

    <input checked=""  disabled="" type="checkbox">

.. |uncheck| raw:: html

    <input type="checkbox">

.. |uncheck_| raw:: html

    <input disabled="" type="checkbox">

Run ventilation code
======================
We use imaging-based biophysical airway tree and ventilation model to simulate airflow in human airways.


Requirements
------------

|uncheck| Python 3.6

|uncheck| Virtual enviornment

|uncheck| Source folder **Packages** - **functional-model**.

|uncheck| CMGUI: a tool for visualising. (**Version 2.7** is highly recommanded by the prime.)

After installation, add the following line to your **.bashrc** file::


   alias cmgui2p7='/hpc/<your upi>/Programs/cmgui-wx-2.7.0_x86_64-linux/cmgui-wx'

Run the ventilation code
------------------------------------------

Step 1: Go to **functional-models** - **ventilation_Swan2011** folder::

   cd functional-models/Ventilation_Swan2011

`FAQ: I am not sure if we need 'grown_full.ipnode and ipelem' to run?`

Step 2: Launch virtual environment


   (Please read the instruction *Python3: Creating a virtual enviornment* if you have not exposed to it.)

Step 3: Run **ventilation_Swan2011.py**

::
 
   python ventilation_Swan2011.py

Step 4: Go to **output** folder. View ventilation output with following:


::

   cd output
   <your upi>@<your computer>:/hpc/<your upi>/Packages/functional-models/ventilation_Swan2011/output$ cmgui2p7

Step 5: In the **Cmgui Command Window**, please try the following:

 1.  Click on: **File** - **Read** - **Node File** - **<P2BRP268-H12816_terminal.exnode>**

 2. Click on: **Graphics** - **Scene Editor**.
	 Then: at ``lines``, choose ``node_points``, click ``add``.
	       at ``Glyph``, choose ``sphere``
	       at ``Base glyph size``, type ``3*3*3`` instead ``1*1*1``.
 3. Back to **Cmgui Command Window**, 
     Click on: **Graphics** - **Spectrum Editor**.
     Click on: **autorange** - **OK**.
 

(Alternatively, if you are confident with the **.com** file. Open it either in the terminal/cmgui, and run.)


