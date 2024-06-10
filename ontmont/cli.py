import click

@click.command()
@click.option('-i', '--input', required=True, help="input file")
def cli(input):
    print(input)

if __name__ == "__main__":
    cli()    
