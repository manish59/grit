class Grit < Formula
  desc "GRIT: Genomic Range Interval Toolkit for genomic interval operations"
  homepage "https://github.com/manish59/grit"
  url "https://github.com/manish59/grit/archive/refs/tags/v0.1.0.tar.gz"
  sha256 "PLACEHOLDER_SHA256"
  license "MIT"
  head "https://github.com/manish59/grit.git", branch: "main"

  depends_on "rust" => :build

  def install
    system "cargo", "install", *std_cargo_args
  end

  test do
    assert_match "grit", shell_output("#{bin}/grit --version")
  end
end
