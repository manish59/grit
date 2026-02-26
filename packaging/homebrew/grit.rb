class Grit < Formula
  desc "GRIT: Genomic Range Interval Toolkit - high-performance BED file operations"
  homepage "https://github.com/manish59/grit"
  url "https://github.com/manish59/grit/archive/refs/tags/v0.1.1.tar.gz"
  sha256 "40c0ad8652a05056aa8eac7662c23d68e67a86343ef588b397f46221c840121d"
  license "MIT"
  head "https://github.com/manish59/grit.git", branch: "main"

  depends_on "rust" => :build

  def install
    system "cargo", "install", *std_cargo_args
  end

  test do
    # Create test BED file
    (testpath/"test.bed").write <<~EOS
      chr1\t100\t200
      chr1\t150\t250
    EOS

    # Test merge command
    output = shell_output("#{bin}/grit merge -i #{testpath}/test.bed")
    assert_match "chr1\t100\t250", output

    # Test version
    assert_match "grit", shell_output("#{bin}/grit --version")
  end
end
