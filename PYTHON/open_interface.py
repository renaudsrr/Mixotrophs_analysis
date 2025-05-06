import tkinter as tk
from tkinter import ttk
import tkinter.font as tkFont
import sys

def open_interface():
    param_dict = {}

    dataset_map = {"LP-NLA": "LN", "CARBBAS": "C", "ELA": "ELA"}
    div_map = {
        "Genus richness": "rich_genus",
        "Shannon diversity": "shannon",
        "Simpson diversity": "op_simpson",
        "Pielou evenness": "eveness_piel",
        "Local Contrib. Beta Div" : "log_LCBD"
    }

    def on_close():
        root.destroy()
        sys.exit(0) 

    def update_test_size_label(value):
        test_size_label.config(text=f"{float(value):.2f}")

    def make_slider_with_label(label_text, variable, from_, to_, is_int=False, help_text=""):
        label = ttk.Label(root, text=label_text)
        slider = ttk.Scale(root, from_=from_, to=to_, variable=variable, orient="horizontal")
        val_label = ttk.Label(root, text=str(variable.get()) if is_int else f"{variable.get():.3f}")
        var_entry = ttk.Entry(root, width=6)
        var_entry.insert(0, str(variable.get()))

        def on_slider_move(v):
            val = int(float(v)) if is_int else float(v)
            val_label.config(text=str(val) if is_int else f"{val:.3f}")
            var_entry.delete(0, tk.END)
            var_entry.insert(0, str(val))

        def on_entry_change(event):
            try:
                val = float(var_entry.get())
                if is_int: val = int(val)
                if from_ <= val <= to_:
                    variable.set(val)
                    slider.set(val)
            except ValueError:
                pass

        slider.config(command=on_slider_move)
        var_entry.bind("<Return>", on_entry_change)
        var_entry.bind("<FocusOut>", on_entry_change)
        help_label = ttk.Label(root, text=help_text, wraplength=300, foreground="gray")
        return (label, slider, val_label, var_entry, help_label)

    def launch():
        param_dict.update({
            "choix_dataset": dataset_map[dataset_display_var.get()],
            "choix_div": div_map[div_display_var.get()],
            "test_size": round(float(test_size_var.get()), 2),
            "harmonisation_dataset": "Y" if harmonisation_var.get() else "N",
            "learning_rate": round(float(learning_rate_var.get()), 4),
            "max_depth": int(max_depth_var.get()),
            "n_estimators": int(n_estimators_var.get()),
            "colsample_bytree": round(float(colsample_bytree_var.get()), 2),
            "subsample": round(float(subsample_var.get()), 2),
            "reg_alpha": round(float(reg_alpha_var.get()), 2),
            "reg_lambda": round(float(reg_lambda_var.get()), 2),
            "min_child_weight": int(min_child_weight_var.get()),
            "gamma": round(float(gamma_var.get()), 2),
            "run_gam_model": "Y" if run_gam_var.get() else "N",
            "correction_div" : "Y" if correction_div.get() else "N",
        })
        root.destroy()

    root = tk.Tk()
    root.title("Model configuration")
    root.option_add("*Font", tkFont.Font(family="Helvetica", size=15))
    root.protocol("WM_DELETE_WINDOW", on_close)

    dataset_display_var = tk.StringVar(value="LP-NLA")
    div_display_var = tk.StringVar(value="Shannon diversity")
    test_size_var = tk.DoubleVar(value=0.2)
    harmonisation_var = tk.BooleanVar(value=True)
    run_gam_var = tk.BooleanVar(value=False)
    correction_div = tk.BooleanVar(value=True)  

    learning_rate_var = tk.DoubleVar(value=0.01)
    max_depth_var = tk.IntVar(value=4)
    n_estimators_var = tk.IntVar(value=1000)
    colsample_bytree_var = tk.DoubleVar(value=0.8)
    subsample_var = tk.DoubleVar(value=0.8)
    reg_alpha_var = tk.DoubleVar(value=0.5)
    reg_lambda_var = tk.DoubleVar(value=3.0)
    min_child_weight_var = tk.IntVar(value=10)
    gamma_var = tk.DoubleVar(value=0.3)

    ttk.Label(root, text="Dataset").grid(row=0, column=0, sticky="w")
    ttk.Combobox(root, textvariable=dataset_display_var, values=list(dataset_map.keys())).grid(row=0, column=1, sticky="w")

    ttk.Label(root, text="Diversity index").grid(row=1, column=0, sticky="w")
    ttk.Combobox(root, textvariable=div_display_var, values=list(div_map.keys())).grid(row=1, column=1, sticky="w")

    ttk.Label(root, text="Test size").grid(row=2, column=0, sticky="w")

    test_slider = ttk.Scale(root, variable=test_size_var, from_=0.05, to=0.5, orient="horizontal", length=150)
    test_slider.config(command=update_test_size_label)
    test_slider.grid(row=2, column=1, sticky="ew")

    test_size_label = ttk.Label(root, text=f"{test_size_var.get():.2f}")
    test_size_label.grid(row=2, column=2, sticky="w", padx=10)

    run_gam_check = ttk.Checkbutton(root, text="Re-RUN GAM", variable=run_gam_var)
    run_gam_check.grid(row=3, column=0, columnspan=2, sticky="w")

    ttk.Checkbutton(root, text="Apply data harmonization", variable=harmonisation_var).grid(row=4, column=0, columnspan=2, sticky="w")
    ttk.Checkbutton(root, text="Apply diversity Autotrophs", variable=correction_div).grid(row=5, column=0, columnspan=2, sticky="w")
    # ttk.Checkbutton(root, text="Apply diversity Autotrophs (-cyano)", variable=correction_div).grid(row=5, column=0, columnspan=2, sticky="w")   


    slider_widgets = [
        make_slider_with_label("Learning rate", learning_rate_var, 0.001, 1.0, help_text="Vitesse d'apprentissage. Plus bas = plus lent mais plus stable"),
        make_slider_with_label("Max depth", max_depth_var, 1, 15, is_int=True, help_text="Profondeur max des arbres (complexité)"),
        make_slider_with_label("n_estimators", n_estimators_var, 100, 3000, is_int=True, help_text="Nombre total d'arbres à entraîner"),
        make_slider_with_label("Colsample_bytree", colsample_bytree_var, 0.0, 1.0, help_text="Proportion de colonnes utilisées par arbre"),
        make_slider_with_label("Subsample", subsample_var, 0.0, 1.0, help_text="Proportion d’échantillons utilisée par arbre"),
        make_slider_with_label("Reg_alpha", reg_alpha_var, 0.0, 5.0, help_text="Régularisation L1 (effet Lasso, pousse à la parcimonie)"),
        make_slider_with_label("Reg_lambda", reg_lambda_var, 0.0, 5.0, help_text="Régularisation L2 (effet Ridge, empêche les gros coeffs)"),
        make_slider_with_label("Min_child_weight", min_child_weight_var, 1, 20, is_int=True, help_text="Poids minimal dans un nœud pour autoriser un split"),
        make_slider_with_label("Gamma", gamma_var, 0.0, 1.0, help_text="Seuil de gain minimal pour autoriser une division de nœud")
    ]

    row = 6
    for widgets in slider_widgets:
        for col, w in enumerate(widgets):
            w.grid(row=row, column=col, sticky="w", padx=4)
        row += 1

    ttk.Button(root, text="Launch", command=launch).grid(row=row, column=0, columnspan=5, pady=10, sticky="ew")
    root.mainloop()
    return param_dict