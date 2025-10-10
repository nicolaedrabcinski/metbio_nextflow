process INSTALL_MINIMAP2 {
    tag "Installing minimap2"
    
    output:
    path "minimap2", emit: minimap2_bin
    
    script:
    """
    echo "📥 Downloading minimap2..."
    
    # Download minimap2 binary
    wget -q https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
    tar -xjf minimap2-2.24_x64-linux.tar.bz2
    
    # Move binary to expected location
    mv minimap2-2.24_x64-linux/minimap2 ./minimap2
    chmod +x minimap2
    
    # Test installation
    echo "🧪 Testing minimap2 installation..."
    if ./minimap2 --version; then
        echo "✅ minimap2 installation successful!"
    else
        echo "❌ ERROR: minimap2 test failed!"
        exit 1
    fi
    """
}