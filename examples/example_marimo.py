import marimo

__generated_with = "0.10.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md("""# Current-carrying capacity calculation for DC cables""")
    return


@app.cell
def _():
    import marimo as mo
    from c4_open.cables import Cable
    from c4_open.project import Project
    from c4_open.arrangements import create_arrangement
    import c4_open
    import pandas as pd
    import matplotlib.pyplot as plt
    return Cable, Project, c4_open, create_arrangement, mo, pd, plt


@app.cell
def _(mo):
    mo.md("""## Project data""")
    return


@app.cell
def _(mo):
    name = mo.ui.text(value='Example', label='Project name')
    name
    return (name,)


@app.cell
def _(mo):
    calculation_case = mo.ui.dropdown(
        options=["dc_sc", "ac_sc", "ac_mc"], value="dc_sc", label='Calculation case')
    calculation_case
    return (calculation_case,)


@app.cell
def _(c4_open):
    cable_type_list = c4_open.cables.get_cable_types()
    return (cable_type_list,)


@app.cell
def _(cable_type_list, mo):
    cable_type = mo.ui.dropdown(
        options=cable_type_list,
        value='A2XH 1x240', label='Cable type')
    cable_type
    return (cable_type,)


@app.cell
def _(mo):
    system_arrangement = mo.ui.dropdown(
        options=["two_flat", "three_flat", "single", "square", "trefoil"], value="two_flat", label='System arrangement')
    system_arrangement
    return (system_arrangement,)


@app.cell
def _(mo):
    L = mo.ui.slider(0.1, 5, 0.1, value=1, label="Depth of laying")
    #N = mo.ui.slider(1, 20, 1, value=2, label='Number of cables')
    number_of_systems = mo.ui.slider(1, 20, 1, value=2, label='Number of systems')
    d_clear = mo.ui.slider(0.01, 1, 0.01, value=0.07, label='Clear distance between systems')
    mo.hstack([L, number_of_systems, d_clear], justify='start')
    return L, d_clear, number_of_systems


@app.cell
def _(mo):
    deltatheta_x = mo.ui.slider(15, 40, 1, value=15, label='Critical temperature rise of soil')
    rho_T4 = mo.ui.slider(0.5, 2.5, 0.5, value=1.0, label='Soil thermal resistivity')
    theta_amb = mo.ui.slider(5, 30, 1, value=20, label='Ambient temperature')
    mo.hstack([deltatheta_x, rho_T4, theta_amb], justify='start')
    return deltatheta_x, rho_T4, theta_amb


@app.cell
def _():
    number_of_cables_per_system = {"two_flat": 2, "three_flat": 3, "single": 1, "square": 4, "trefoil": 3}
    return (number_of_cables_per_system,)


@app.cell
def _(number_of_cables_per_system, number_of_systems, system_arrangement):
    N = number_of_systems.value * number_of_cables_per_system[system_arrangement.value]
    return (N,)


@app.cell
def _(
    L,
    N,
    Project,
    cable_type,
    calculation_case,
    create_arrangement,
    d_clear,
    deltatheta_x,
    name,
    number_of_systems,
    rho_T4,
    system_arrangement,
    theta_amb,
):
    project_data = {'name': name.value, 
                    'calc_case': calculation_case.value, 
                    'cable_type':cable_type.value,
                    'L': L.value,
                    'N': N,
                    'deltatheta_x': deltatheta_x.value,
                    'rho_T4': rho_T4.value,
                    'theta_amb': theta_amb.value
                   }

    project = Project(**project_data)

    F,_,xx,yy = create_arrangement(system_arrangement.value, number_of_systems=number_of_systems.value, d_ext=project.D_e / 1000, d_clear=d_clear.value, depth_of_laying=project.L)

    project.F = F
    return F, project, project_data, xx, yy


@app.cell
def _(plt, xx, yy):
    plt.scatter(xx, yy, marker='o', s = 3)
    plt.xlim(-1, 1)
    plt.ylim(-1.5, 0.7)
    plt.gca()
    return


@app.cell
def _(Cable, project):
    cable = Cable(project)
    return (cable,)


@app.cell
def _(cable):
    result = cable.get_result()
    return (result,)


@app.cell
def _(mo):
    mo.md("""## Resulting current carrying capacity""")
    return


@app.cell
def _(mo, result):
    mo.output.replace(f'{round(result['I (with dry zone)'], 1)} A')
    return


if __name__ == "__main__":
    app.run()
