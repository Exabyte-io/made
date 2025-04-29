import os
import tempfile
import time
import webbrowser

# Use a default div id
default_div_id = "wave-div"


def get_wave_html(div_id=default_div_id, width=600, height=600, title="Material"):
    size = min(width, height)  # Make it square using the smaller dimension
    return f"""
    <h2>{title}</h2>
    <div id="{div_id}" style="width:{size}px; height:{size}px; border:1px solid #333;"></div>
    """


def get_wave_js(material_json, div_id=default_div_id):
    return (
        f"""
    const materialConfig = {material_json};
    const container = document.getElementById('{div_id}');
        """
        + """
    (async function() {
        const module = await import('https://exabyte-io.github.io/wave.js/main.js');
        window.renderThreeDEditor(materialConfig, container);
    })();
    document.head.insertAdjacentHTML(
        'beforeend',
        '<link rel="stylesheet" href="https://exabyte-io.github.io/wave.js/main.css"/>'
    );
    """
    )


def debug_visualize_material(material, width=600, height=600, title="Material"):
    """
    Generates a temporary HTML file that uses Wave.js to visualize the material,
    and opens it in the default browser.

    Call this function from the PyCharm debugger (e.g., via Evaluate Expression).
    """
    # Convert your material to JSON.
    # (Assuming material.to_json() returns a JSON-serializable object)
    material_json = material.to_json()

    # Generate a unique div id so multiple calls don't conflict
    div_id = f"wave-{int(time.time())}"

    # Create HTML content that includes our working code
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="UTF-8">
      <title>Wave.js Debug Viewer</title>
    </head>
    <body>
      {get_wave_html(div_id, width, height, title)}
      <script type="module">
        {get_wave_js(material_json, div_id)}
      </script>
    </body>
    </html>
    """

    # Write the HTML to a temporary file and open it in the default browser
    fd, file_path = tempfile.mkstemp(suffix=".html", prefix="wave_debug_")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write(html_content)

    webbrowser.open("file://" + file_path)
