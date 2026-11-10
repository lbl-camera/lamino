import click
import numpy as np
import matplotlib.pyplot as plt


def plot_(csvfile, output='convergence_plot.png'):

    # get the header from the csv file
    with open(csvfile, 'r') as f:
        header = f.readline().strip().split(',')

    # load the data from the csv file, skipping the header
    data = np.loadtxt(csvfile, delimiter=',', skiprows=1)

    # create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(data[:, 0], data[:, 1], label=header[1], color='blue', marker='o')
    ax.set_xlabel(header[0])
    ax.set_ylabel('Values')
    # right side y-axis for the second column
    ax2 = ax.twinx()
    ax2.semilogy(data[:, 0], data[:, 2], label=header[2], color='green', marker='o')
    ax2.set_ylabel(header[2])
    ax.set_title('Convergence Plot')
    ax.grid()

    # add a legend
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc='best')

    plt.tight_layout()
    plt.savefig(output)
    plt.show()


@click.command()
@click.argument('csvfile', type=click.Path(exists=True))
@click.option('--output', '-o', default='convergence_plot.png', help='Output file name for the plot')
def plot_convergence(csvfile, output):
    """Plot convergence from a CSV file."""
    plot_(csvfile, output)
    click.echo(f'Convergence plot saved as {output}')


if __name__ == '__main__':
    plot_convergence()
