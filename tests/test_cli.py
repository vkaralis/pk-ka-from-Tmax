from ka_from_tmax_kel.cli import main


def test_cli_output(capsys):
    assert main(["1.5", "0.5"]) == 0
    output = capsys.readouterr().out
    assert "Ka = 0.86680049" in output
    assert "Residual Tmax" in output
