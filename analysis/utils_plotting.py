from scipy import stats
from matplotlib.cm import get_cmap
import seaborn as sns


def get_color_palette(num):
    """
    Returns color palette to use for plotting emerging lineages. Colors match those used by Nextstrain.
    """
    color_palettes = {8:['#4068CF','#5098B9','#6CB28C','#94BD62','#BFBB47','#DFA53B','#E67131','#DB2823'],
    9: ['#3E5DD0','#4A8CC2','#60AA9E','#80B974','#A6BE55','#CBB742','#E29D39','#E56A2F','#DB2823'],
    10:['#3F52CD','#4681CA','#57A1AD','#70B487','#90BC65','#B4BD4C','#D3B240','#E59638','#E4642E','#DB2823'],
    11:['#3F47C9','#4274CE','#4F97BB','#64AC99','#7EB976','#9EBE5A','#BEBB48','#D9AE3E','#E69036','#E35F2D','#DB2823'],
    12:['#403CC5','#4067CF','#4A8BC3','#5AA4A9','#6FB488','#8BBB69','#A9BD53','#C7B944','#DDA93C','#E68A35','#E35C2C','#DB2823'],
    13:['#4433BE','#3E5ACF','#457FCB','#529AB6','#64AD98','#7BB77A','#96BD60','#B3BD4D','#CDB642','#DFA43B','#E68434','#E2582C','#DB2823'],
    14:['#492AB5','#3F4CCB','#4271CE,''#4C8FC0','#5AA5A8','#6DB38A','#85BA6F','#A0BE59','#BBBC49','#D2B340','#E19F3A','#E68033','#E2562B','#DB2823'],
    15:['#4D21AD','#403FC6','#3F63CF','#4783C8','#539BB5','#63AC9B','#77B67F','#8EBC66','#A8BD53','#C1BA47','#D6B03F','#E39C39','#E67C33','#E1532B','#DB2823'],
    16:['#571EA2','#4334BF','#3F55CE','#4376CD','#4C91C0','#59A4A9','#6AB18F','#7FB975','#97BD5F','#AFBD4F','#C7B944','#D9AD3D','#E49838','#E67932','#E1512A','#DB2823'],
    17:['#5E1D9D','#462EB9','#3F4CCB','#416CCE','#4887C6','#539CB3','#62AB9C','#74B582','#89BB6B','#A0BE59','#B7BD4B','#CCB742','#DDAA3C','#E69537','#E67631','#E14F2A','#DB2823']}

    return color_palettes[num]

def convert_linege_names(old_name):
    """
    Converts old names of emerging lineages into new names, including WHO VOC designation. Also renames 'unassigned' to 'basal'
    """
    old_to_new = {'A.23.1': 'A.23.1', 'B.1.1.7': 'B.1.1.7 (Alpha)', 'B.1.351': 'B.1.351 (Beta)',
                  'B.1.427+B.1.429': 'B.1.427/429 (Epsilon)', 'B.1.525': 'B.1.525 (Eta)', 'B.1.526': 'B.1.526 (Iota)',
                  'B.1.617': 'B.1.617.1/2 (Delta/Kappa)',
                  'C.37': 'C.37 (Lambda)', 'P.1': 'P.1 (Gamma)', 'P.3': 'P.3',
                  'unassigned': 'basal'}

    if old_name in old_to_new.keys():
        new_name = old_to_new[old_name]
    else:
        new_name = old_name

    return new_name

def hue_regplot(data, x, y, hue, palette=None, **kwargs):
    """
    To plot seaborn regpot with a hue option because lmplot doesn't have 'ax'
    """

    regplots = []

    levels = data[hue].unique()

    if palette is None:
        default_colors = get_cmap('tab10')
        palette = {k: default_colors(i) for i, k in enumerate(levels)}

    for key in levels:
        regplots.append(
            sns.regplot(
                x=x,
                y=y,
                data=data[data[hue] == key],
                color=palette[key],
                **kwargs
            )
        )

    return regplots

def get_linear_reg_stats(df, mut_location, x_axis, y_axis):

    df_subset = df.copy()
    df_subset = df_subset[df_subset['mut_location']==mut_location]

    if x_axis == 'logistic_growth':
        df_subset = df_subset.dropna()


    slope, intercept, r_value, p_value, std_err = stats.linregress(df_subset[x_axis], df_subset[y_axis])
    slope = round(slope, 4)
    r_value = round(r_value, 2)

    return slope, r_value
